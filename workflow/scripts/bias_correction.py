import os
import logging

import numpy as np
import pandas as pd

import statsmodels.formula.api as smf
import statsmodels.api as sm_ols
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
import patsy

from pybedtools import BedTool
import pyranges as pr
import matplotlib.pyplot as plt


def compute_gc_content(
    bin_info: pd.DataFrame, ref_file: str, mapp_file=None, genome_size=None
):
    """Compute per-bin GC fraction and optionally mean mappability.

    Uses pybedtools ``nucleotide_content`` on the reference FASTA. If a
    mappability BED track is provided, per-bin mean mappability is also computed.

    Parameters
    ----------
    bin_info : pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` columns defining genomic bins.
    ref_file : str
        Path to the reference genome FASTA.
    mapp_file : str or None
        Path to a BED-format mappability track (optional).
    genome_size : str or None
        Path to chromosome sizes file (required when *mapp_file* is given).

    Returns
    -------
    pd.DataFrame
        *bin_info* augmented with ``GC`` and ``MAP`` columns.
    """
    logging.info("compute GC content")
    gc_df = bin_info.merge(
        BedTool.from_dataframe(bin_info[["#CHR", "START", "END"]])
        .nucleotide_content(fi=ref_file)
        .to_dataframe(disable_auto_names=True)
        .rename(
            columns={
                "#1_usercol": "#CHR",
                "2_usercol": "START",
                "3_usercol": "END",
                "5_pct_gc": "GC",
            }
        )[["#CHR", "START", "END", "GC"]]
    )
    if mapp_file:
        logging.info("compute mappability")
        map_bt = BedTool(mapp_file)
        map_cov = (
            BedTool.from_dataframe(bin_info[["#CHR", "START", "END"]])
            .map(b=map_bt, c=4, o="mean", g=genome_size)
            .to_dataframe(disable_auto_names=True)
        )
        map_cov.columns = ["#CHR", "START", "END", "MAP"]
        map_cov["MAP"] = (
            pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(1.0).clip(0.0, 1.0)
        )
        gc_df = gc_df.merge(map_cov, on=["#CHR", "START", "END"], how="left")
    else:
        gc_df["MAP"] = 1.0
    return gc_df


def bias_correction_rdr(
    raw_rdr_mat: np.ndarray,
    gc_df: pd.DataFrame,
    rep_ids: list,
    has_mapp=False,
    out_dir=None,
    eps_quantile=0.01,
    gc_quantile=[0.01, 0.99],
):
    """Correct read-depth ratios (RDR) for GC-content (and optionally mappability) bias.

    Fits a median quantile regression of RDR on GC (quadratic) per tumor sample,
    divides raw RDR by the fitted expected RDR, and normalizes by the mean
    corrected RDR. Saves diagnostic scatter and hexbin plots.

    Parameters
    ----------
    raw_rdr_mat : np.ndarray
        Raw RDR matrix (bins x samples).
    gc_df : pd.DataFrame
        DataFrame with ``GC`` and ``MAP`` columns per bin.
    rep_ids : list[str]
        Replicate identifiers (tumor samples only).
    has_mapp : bool
        If True, include mappability as a covariate.
    out_dir : str or None
        Directory for diagnostic plots.
    eps_quantile : float
        Lower quantile of expected RDR used as floor to avoid division by zero.
    gc_quantile : list[float]
        Lower and upper GC quantiles defining the fitting range.

    Returns
    -------
    np.ndarray
        GC-corrected RDR matrix (bins x samples).
    """
    logging.info("correct for GC biases on RDR")
    gc = gc_df["GC"].to_numpy()
    mapv = gc_df["MAP"].to_numpy()

    gc_lo, gc_hi = np.nanquantile(gc, gc_quantile)
    logging.info(f"GC-content range for fitting: [{gc_lo}, {gc_hi}]")
    fit_mask = (gc >= gc_lo) & (gc <= gc_hi)

    gccorr_rdr_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        raw_rdrs = raw_rdr_mat[:, si]
        gc_df = pd.DataFrame({"RD": raw_rdrs[fit_mask], "GC": gc[fit_mask]})
        if has_mapp:
            gc_df["MAP"] = mapv[fit_mask]
            res = smf.quantreg("RD ~ GC + I(GC**2) + MAP + I(MAP**2)", data=gc_df).fit(
                q=0.5
            )
            logging.info(res.summary())
            exp_rdrs = res.predict(pd.DataFrame({"GC": gc, "MAP": mapv}))
        else:
            res = smf.quantreg("RD ~ GC + I(GC**2)", data=gc_df).fit(q=0.5)
            logging.info(res.summary())
            exp_rdrs = res.predict(pd.DataFrame({"GC": gc})).to_numpy()
        eps = np.nanquantile(exp_rdrs, eps_quantile)
        logging.info(f"inferred epsilon for exp RDR={eps}")
        den = exp_rdrs.copy()
        den[~np.isfinite(den)] = np.nan
        den = np.where(np.isfinite(den), den, eps)
        den = np.clip(den, eps, None)
        corr_rdrs = raw_rdrs / den

        corr_factor = np.mean(corr_rdrs)
        logging.info(f"GC correction factor={corr_factor}")

        # raw RDRs
        logging.info("hist: raw RDR")
        counts, edges = np.histogram(raw_rdrs)
        logging.info("bin_left\tbin_right\tcount")
        for l, r, c in zip(edges[:-1], edges[1:], counts):
            logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")

        # expected RDRs
        logging.info("hist: expected RDR corr_fit")
        counts, edges = np.histogram(exp_rdrs)
        logging.info("bin_left\tbin_right\tcount")
        for l, r, c in zip(edges[:-1], edges[1:], counts):
            logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")

        logging.info("hist: corr_rdrs")
        counts, edges = np.histogram(corr_rdrs)
        logging.info("bin_left\tbin_right\tcount")
        for l, r, c in zip(edges[:-1], edges[1:], counts):
            logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")

        gccorr_rdr_mat[:, si] = corr_rdrs / corr_factor
        if out_dir is not None:
            plot_correction_diagnostics(
                gc, raw_rdrs, gccorr_rdr_mat[:, si],
                fitted_rdr=exp_rdrs, rep_id=rep_id, out_dir=out_dir,
            )
    return gccorr_rdr_mat


def plot_correction_diagnostics(
    gc, raw_rdr, corr_rdr, fitted_rdr, rep_id, out_dir,
    rt_vals=None, rt_name=None,
):
    """Produce all diagnostic plots for RDR bias correction.

    Generates three plot files:

    1. ``{rep_id}.gc_corr_scatter.png`` — GC vs raw RDR scatter with fitted
       expected RDR overlaid as a red line.
    2. ``{rep_id}.gc_corr.png`` — side-by-side hexbin of uncorrected vs.
       corrected RDR, each annotated with MAD.
    3. ``{rep_id}.rt_scatter.png`` — (only when *rt_vals* is provided)
       before/after RT correction scatter.

    Parameters
    ----------
    gc : np.ndarray
        Per-bin GC content.
    raw_rdr : np.ndarray
        Uncorrected RDR values.
    corr_rdr : np.ndarray
        Corrected RDR values.
    fitted_rdr : np.ndarray
        Model-predicted expected RDR for every bin (used for the fit overlay).
    rep_id : str
        Replicate identifier (used in output filenames).
    out_dir : str
        Directory for the output PNGs.
    rt_vals : np.ndarray or None
        Per-bin replication timing values (optional).
    rt_name : str or None
        Name of the selected RT cell line (used in axis labels).
    """
    def mad(x):
        return np.median(np.abs(x - np.median(x)))

    # --- 1. GC vs RDR scatter with fit ---
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.scatter(gc, raw_rdr, s=2, alpha=0.2, label="raw RDR")
    sort_idx = np.argsort(gc)
    ax.plot(
        gc[sort_idx], fitted_rdr[sort_idx],
        linewidth=2, color="red", label="fit",
    )
    ax.set_xlabel("GC")
    ax.set_ylabel("RDR")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{rep_id}.gc_corr_scatter.png"), dpi=300)
    plt.close(fig)

    # --- 2. Before/after hexbin ---
    mad_raw = mad(raw_rdr)
    mad_corr = mad(corr_rdr)
    fig, axes = plt.subplots(1, 2, figsize=(6, 3), sharey=True)
    axes[0].hexbin(
        gc, raw_rdr,
        gridsize=100, cmap="Blues", mincnt=1, linewidths=0,
        reduce_C_function=np.median,
    )
    axes[0].set_xlabel("GC content")
    axes[0].set_ylabel("RDR")
    axes[0].set_title(f"Uncorrected RDR\nMAD={mad_raw:.3f}")
    axes[1].hexbin(
        gc, corr_rdr,
        gridsize=100, cmap="Blues", mincnt=1, linewidths=0,
        reduce_C_function=np.median,
    )
    axes[1].set_xlabel("GC content")
    axes[1].set_ylabel("RDR")
    axes[1].set_title(f"GC-corrected RDR\nMAD={mad_corr:.3f}")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{rep_id}.gc_corr.png"), dpi=300)
    plt.close(fig)

    # --- 3. RT scatter (optional) ---
    if rt_vals is not None:
        valid = np.isfinite(rt_vals)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].scatter(rt_vals[valid], raw_rdr[valid], s=2, alpha=0.2)
        axes[0].set_xlabel(f"RT ({rt_name})")
        axes[0].set_ylabel("RDR (uncorrected)")
        axes[0].set_title("Before RT correction")
        axes[1].scatter(rt_vals[valid], corr_rdr[valid], s=2, alpha=0.2)
        axes[1].set_xlabel(f"RT ({rt_name})")
        axes[1].set_ylabel("RDR (corrected)")
        axes[1].set_title("After RT correction")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{rep_id}.rt_scatter.png"), dpi=300)
        plt.close(fig)


def load_rt_for_bins(bin_info: pd.DataFrame, rt_file: str) -> pd.DataFrame:
    """Load ASCAT replication timing data and average over genomic bins.

    Reads a tab-separated ASCAT RT file (columns: Chr, Position, then one column
    per cell line) and maps SNP-level RT values to bins via pyranges overlap join,
    averaging RT within each bin.

    Parameters
    ----------
    bin_info : pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` columns defining genomic bins.
    rt_file : str
        Path to ASCAT replication timing TSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` and one column per RT cell line.
    """
    logging.info(f"loading replication timing from {rt_file}")
    rt_raw = pd.read_csv(rt_file, sep="\t", index_col=0)
    # Columns: Chr, Position, <cell_line_1>, ..., <cell_line_N>
    rt_raw["Chr"] = rt_raw["Chr"].astype(str).str.replace("^chr", "", regex=True)
    cell_line_cols = [c for c in rt_raw.columns if c not in ("Chr", "Position")]

    # Build pyranges for RT SNPs (point intervals)
    rt_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": rt_raw["Chr"],
        "Start": rt_raw["Position"].astype(int),
        "End": rt_raw["Position"].astype(int) + 1,
        **{c: rt_raw[c].values for c in cell_line_cols},
    }))

    # Build pyranges for bins
    bin_pr = pr.PyRanges(pd.DataFrame({
        "Chromosome": bin_info["#CHR"].astype(str),
        "Start": bin_info["START"].astype(int),
        "End": bin_info["END"].astype(int),
    }))

    # Overlap join and aggregate by bin
    joined = bin_pr.join(rt_pr, how="left")
    joined_df = joined.df if hasattr(joined, "df") else joined.as_df()

    # Group by bin coordinates, take mean of each cell line
    agg = (
        joined_df
        .groupby(["Chromosome", "Start", "End"])[cell_line_cols]
        .mean()
        .reset_index()
    )

    # Merge back to preserve bin order
    result = (
        bin_info[["#CHR", "START", "END"]]
        .merge(
            agg.rename(columns={"Chromosome": "#CHR", "Start": "START", "End": "END"}),
            on=["#CHR", "START", "END"],
            how="left",
        )
    )

    n_missing = result[cell_line_cols[0]].isna().sum()
    n_total = len(result)
    frac_missing = n_missing / n_total
    if frac_missing > 0.05:
        logging.warning(
            f"RT: {n_missing}/{n_total} bins ({frac_missing:.1%}) have no overlapping "
            f"RT SNPs; these will be filled with NaN and excluded from fitting."
        )
    logging.info(f"RT loaded for {len(cell_line_cols)} cell lines, {n_total} bins")
    return result


def select_best_rt_cellline(
    rdr_vec: np.ndarray, rt_df: pd.DataFrame, autosomes_mask: np.ndarray
) -> tuple:
    """Select the RT cell line most correlated with RDR on autosomal bins.

    Parameters
    ----------
    rdr_vec : np.ndarray
        Per-bin RDR values for one sample.
    rt_df : pd.DataFrame
        DataFrame with RT cell line columns (from ``load_rt_for_bins``).
    autosomes_mask : np.ndarray
        Boolean mask for autosomal bins (excludes X, Y).

    Returns
    -------
    tuple[str, float]
        (best_cellline_name, absolute_pearson_correlation).
    """
    cell_line_cols = [c for c in rt_df.columns if c not in ("#CHR", "START", "END")]
    best_name, best_corr = None, -1.0

    for cl in cell_line_cols:
        rt_vals = rt_df[cl].to_numpy()
        valid = autosomes_mask & np.isfinite(rdr_vec) & np.isfinite(rt_vals)
        if valid.sum() < 10:
            logging.info(f"RT cell line {cl}: too few valid bins ({valid.sum()}), skipping")
            continue
        r, _ = pearsonr(rdr_vec[valid], rt_vals[valid])
        abs_r = abs(r)
        logging.info(f"RT cell line {cl}: |r| = {abs_r:.4f}")
        if abs_r > best_corr:
            best_corr = abs_r
            best_name = cl

    logging.info(f"best RT cell line: {best_name} (|r| = {best_corr:.4f})")
    return best_name, best_corr


def bias_correction_rdr_spline(
    raw_rdr_mat: np.ndarray,
    gc_df: pd.DataFrame,
    rep_ids: list,
    has_mapp=False,
    rt_df=None,
    out_dir=None,
    gc_quantile=(0.01, 0.99),
):
    """ASCAT-style spline RDR correction using natural splines on GC (+MAP, +RT).

    Works in log-space: fits OLS on natural cubic spline basis (df=5) of GC and
    optionally mappability and replication timing. Corrected RDR is computed as
    residuals shifted to the mean log-RDR, then exponentiated.

    Parameters
    ----------
    raw_rdr_mat : np.ndarray
        Raw RDR matrix (bins x samples).
    gc_df : pd.DataFrame
        DataFrame with ``GC``, ``MAP``, ``#CHR`` columns per bin.
    rep_ids : list[str]
        Replicate identifiers (tumor samples only).
    has_mapp : bool
        If True, include mappability as a spline covariate.
    rt_df : pd.DataFrame or None
        RT DataFrame from ``load_rt_for_bins``, or None to skip RT.
    out_dir : str or None
        Directory for diagnostic plots.
    gc_quantile : tuple[float, float]
        Lower and upper GC quantiles defining the fitting range.

    Returns
    -------
    np.ndarray
        Corrected RDR matrix (bins x samples).
    """
    logging.info("ASCAT-style spline RDR correction")
    gc = gc_df["GC"].to_numpy()
    mapv = gc_df["MAP"].to_numpy()
    chrs = gc_df["#CHR"].astype(str).to_numpy()

    gc_lo, gc_hi = np.nanquantile(gc, gc_quantile)
    logging.info(f"GC-content range for fitting: [{gc_lo}, {gc_hi}]")

    autosomes_mask = np.array([c not in ("X", "Y", "chrX", "chrY") for c in chrs])

    corrected_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        raw_rdrs = raw_rdr_mat[:, si]
        log_rdr = np.log(raw_rdrs + 1e-8)

        # Select best RT cell line if RT data available
        rt_vals = None
        best_rt_name = None
        if rt_df is not None:
            best_rt_name, best_rt_corr = select_best_rt_cellline(
                raw_rdrs, rt_df, autosomes_mask
            )
            if best_rt_name is not None:
                rt_vals = rt_df[best_rt_name].to_numpy()

        # Build covariate DataFrame
        cov_df = pd.DataFrame({"log_rdr": log_rdr, "GC": gc})
        formula_parts = ["cr(GC, df=5)"]
        if has_mapp:
            cov_df["MAP"] = mapv
            formula_parts.append("cr(MAP, df=5)")
        if rt_vals is not None:
            cov_df["RT"] = rt_vals
            formula_parts.append("cr(RT, df=5)")

        formula_rhs = " + ".join(formula_parts)

        # Determine fit mask: within GC quantile range, finite values
        fit_mask = (
            (gc >= gc_lo) & (gc <= gc_hi) &
            np.isfinite(log_rdr)
        )
        if has_mapp:
            fit_mask &= np.isfinite(mapv)
        if rt_vals is not None:
            fit_mask &= np.isfinite(rt_vals)

        logging.info(f"{rep_id}: fitting spline on {fit_mask.sum()}/{len(fit_mask)} bins")

        # Build design matrix on fit subset
        fit_df = cov_df.loc[fit_mask].copy()
        y_fit = fit_df["log_rdr"].values
        X_fit = patsy.dmatrix(f"{formula_rhs} - 1", data=fit_df, return_type="dataframe")

        # Fit OLS
        ols_model = sm_ols.OLS(y_fit, X_fit).fit()
        logging.info(f"{rep_id}: spline OLS R²={ols_model.rsquared:.4f}")

        # Predict on ALL bins
        # For bins with NaN covariates, fill temporarily for prediction, then mark as NaN
        pred_df = cov_df.copy()
        for col in pred_df.columns:
            if col == "log_rdr":
                continue
            pred_df[col] = pred_df[col].fillna(pred_df[col].median())
        X_all = patsy.dmatrix(
            f"{formula_rhs} - 1", data=pred_df, return_type="dataframe"
        )
        predicted = ols_model.predict(X_all)

        # Residuals
        residuals = log_rdr - predicted
        mean_log_rdr_fit = np.mean(log_rdr[fit_mask])
        corrected_log_rdr = residuals + mean_log_rdr_fit

        # Back-transform
        corrected_rdr = np.exp(corrected_log_rdr)

        # Handle NaN bins (those with NaN in original covariates): impute with median
        nan_mask = ~np.isfinite(log_rdr)
        if rt_vals is not None:
            nan_mask |= ~np.isfinite(rt_vals)
        if nan_mask.any():
            med_corr = np.nanmedian(corrected_rdr[~nan_mask])
            corrected_rdr[nan_mask] = med_corr
            logging.info(
                f"{rep_id}: imputed {nan_mask.sum()} NaN bins with median={med_corr:.4f}"
            )

        # Normalize by per-sample mean
        corr_factor = np.mean(corrected_rdr)
        logging.info(f"{rep_id}: spline correction factor={corr_factor:.4f}")
        corrected_mat[:, si] = corrected_rdr / corr_factor

        # Diagnostic plots
        if out_dir is not None:
            fitted_rdr = np.exp(predicted)
            plot_correction_diagnostics(
                gc, raw_rdrs, corrected_mat[:, si],
                fitted_rdr=fitted_rdr, rep_id=rep_id, out_dir=out_dir,
                rt_vals=rt_vals, rt_name=best_rt_name,
            )

    return corrected_mat
