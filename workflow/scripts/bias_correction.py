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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def _make_autosomes_mask(chrs):
    """Boolean mask that is True for autosomal chromosomes."""
    return np.array([c not in ("X", "Y", "chrX", "chrY") for c in chrs])


def _select_rt_for_sample(raw_rdrs, rt_df, autosomes_mask):
    """Pick best RT cell line for one sample; return (rt_vals, name) or (None, None)."""
    if rt_df is None:
        return None, None
    name, _ = select_best_rt_cellline(raw_rdrs, rt_df, autosomes_mask)
    if name is None:
        return None, None
    return rt_df[name].to_numpy(), name


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


def bias_correction_rdr_quantreg(
    raw_rdr_mat: np.ndarray,
    gc_df: pd.DataFrame,
    rep_ids: list,
    has_mapp=False,
    rt_df=None,
    out_dir=None,
    gc_quantile=[0.05, 0.95],
    eps=1e-8,
):
    """Correct RDR for GC (+ mappability, + RT) bias via median quantile regression in log2 space.

    Fits a quadratic median quantile regression of log2(RDR) on GC per sample,
    with optional mappability and replication-timing covariates. Corrected RDR
    is computed from residuals shifted to the sample mean, then exponentiated
    and mean-normalized to 1.

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
    rt_df : pd.DataFrame or None
        RT DataFrame from ``load_rt_for_bins``, or None to skip RT.
    out_dir : str or None
        Directory for diagnostic plots.
    gc_quantile : list[float]
        Lower and upper GC quantiles defining the fitting range.
    eps : float
        Small constant added before log2 to avoid log(0).

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
    autosomes_mask = _make_autosomes_mask(gc_df["#CHR"].astype(str).to_numpy())

    pdf = PdfPages(os.path.join(out_dir, "correction_diagnostics.pdf")) if out_dir else None

    corrected_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        raw_rdrs = raw_rdr_mat[:, si]
        log_rdr = np.log2(raw_rdrs + eps)
        rt_vals, best_rt_name = _select_rt_for_sample(raw_rdrs, rt_df, autosomes_mask)

        # Build formula, fit_df, pred_df dynamically
        sample_fit_mask = fit_mask.copy()
        sample_fit_mask &= np.isfinite(log_rdr)
        if rt_vals is not None:
            sample_fit_mask &= np.isfinite(rt_vals)

        fit_df = pd.DataFrame({"log_RD": log_rdr[sample_fit_mask], "GC": gc[sample_fit_mask]})
        pred_df = pd.DataFrame({"GC": gc})
        formula = "log_RD ~ GC + I(GC**2)"
        if has_mapp:
            fit_df["MAP"] = mapv[sample_fit_mask]
            pred_df["MAP"] = mapv
            formula += " + MAP + I(MAP**2)"
        if rt_vals is not None:
            fit_df["RT"] = rt_vals[sample_fit_mask]
            pred_df["RT"] = np.where(np.isfinite(rt_vals), rt_vals, np.nanmedian(rt_vals))
            formula += " + RT + I(RT**2)"

        res = smf.quantreg(formula, data=fit_df).fit(q=0.5)
        logging.info(res.summary())
        predicted = res.predict(pred_df).to_numpy()

        residuals = log_rdr - predicted
        mean_log_rdr_fit = np.mean(log_rdr[sample_fit_mask])
        corrected_log_rdr = residuals + mean_log_rdr_fit
        corr_rdrs = np.exp2(corrected_log_rdr)

        nan_mask = ~np.isfinite(log_rdr)
        if rt_vals is not None:
            nan_mask |= ~np.isfinite(rt_vals)
        if nan_mask.any():
            med_corr = np.nanmedian(corr_rdrs[~nan_mask])
            corr_rdrs[nan_mask] = med_corr
            logging.info(
                f"{rep_id}: imputed {nan_mask.sum()} NaN bins with median={med_corr:.4f}"
            )

        exp_rdrs = np.exp2(predicted)
        corrected_mat[:, si] = corr_rdrs / np.mean(corr_rdrs)
        if pdf is not None:
            plot_correction_diagnostics(
                gc, raw_rdrs, corrected_mat[:, si],
                fitted_rdr=exp_rdrs, rep_id=rep_id, pdf=pdf,
                rt_vals=rt_vals, rt_name=best_rt_name,
            )

    if pdf is not None:
        pdf.close()
    return corrected_mat


def plot_correction_diagnostics(
    gc, raw_rdr, corr_rdr, fitted_rdr, rep_id, pdf,
    rt_vals=None, rt_name=None,
):
    """Produce a single-page diagnostic figure for RDR bias correction.

    Adds one page to *pdf* with a 2x2 scatter grid (GC before/after on row 0,
    RT before/after on row 1). If no RT data is provided, only the top row
    (1x2) is drawn.

    Parameters
    ----------
    gc : np.ndarray
        Per-bin GC content.
    raw_rdr : np.ndarray
        Uncorrected RDR values.
    corr_rdr : np.ndarray
        Corrected RDR values.
    fitted_rdr : np.ndarray
        Model-predicted expected RDR for every bin.
    rep_id : str
        Replicate identifier (used in page title).
    pdf : PdfPages
        Open PdfPages handle to write the figure into.
    rt_vals : np.ndarray or None
        Per-bin replication timing values (optional).
    rt_name : str or None
        Name of the selected RT cell line (used in axis labels).
    """
    mad = lambda x: np.median(np.abs(x - np.median(x)))
    fitted_rdr = np.asarray(fitted_rdr)

    has_rt = rt_vals is not None
    if has_rt:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes = axes[np.newaxis, :]  # make 2D for uniform indexing

    # Row 0: GC scatter
    mad_raw = mad(raw_rdr)
    mad_corr = mad(corr_rdr)
    axes[0, 0].scatter(gc, raw_rdr, s=2, alpha=0.2, rasterized=True)
    sort_idx = np.argsort(gc)
    axes[0, 0].plot(gc[sort_idx], fitted_rdr[sort_idx], linewidth=2, color="red", label="fit")
    axes[0, 0].legend()
    axes[0, 0].set_xlabel("GC")
    axes[0, 0].set_ylabel("RDR")
    axes[0, 0].set_title(f"Before GC correction\nMAD={mad_raw:.4f}")
    axes[0, 1].scatter(gc, corr_rdr, s=2, alpha=0.2, rasterized=True)
    axes[0, 1].set_xlabel("GC")
    axes[0, 1].set_ylabel("RDR")
    axes[0, 1].set_title(f"After GC correction\nMAD={mad_corr:.4f}")

    # Row 1: RT scatter (only if RT data available)
    if has_rt:
        valid = np.isfinite(rt_vals)
        mad_raw_rt = mad(raw_rdr[valid])
        mad_corr_rt = mad(corr_rdr[valid])
        axes[1, 0].scatter(rt_vals[valid], raw_rdr[valid], s=2, alpha=0.2, rasterized=True)
        axes[1, 0].scatter(rt_vals[valid], fitted_rdr[valid], s=2, alpha=0.2, color="red", label="fit", rasterized=True)
        axes[1, 0].legend()
        axes[1, 0].set_xlabel(f"RT ({rt_name})")
        axes[1, 0].set_ylabel("RDR")
        axes[1, 0].set_title(f"Before RT correction\nMAD={mad_raw_rt:.4f}")
        axes[1, 1].scatter(rt_vals[valid], corr_rdr[valid], s=2, alpha=0.2, rasterized=True)
        axes[1, 1].set_xlabel(f"RT ({rt_name})")
        axes[1, 1].set_ylabel("RDR")
        axes[1, 1].set_title(f"After RT correction\nMAD={mad_corr_rt:.4f}")

    fig.suptitle(rep_id, fontsize=14)
    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def load_rt_for_bins(bin_info: pd.DataFrame, rt_file: str) -> pd.DataFrame:
    """Load ASCAT replication timing data and average over genomic bins.

    Reads a tab-separated ASCAT RT file (columns: Chr, Position, then one column
    per cell line) and assigns each RT SNP to a bin via ``np.searchsorted``,
    then averages RT values within each bin.

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
    rt_raw = pd.read_csv(rt_file, sep="\t", index_col=None, dtype={"Chr": str, "Position": np.int64})
    rt_raw = rt_raw.drop(columns=[rt_raw.columns[0]])
    logging.info(f"RT file loaded: {len(rt_raw)} SNPs, {rt_raw.shape[1]} columns")
    if not str.startswith(rt_raw["Chr"].iloc[0], "chr"):
        rt_raw["Chr"] = "chr" + rt_raw["Chr"]
    cell_line_cols = [c for c in rt_raw.columns if c not in ("Chr", "Position")]
    logging.info(f"RT cell lines: {cell_line_cols}")

    bin_chr = bin_info["#CHR"].astype(str).values
    bin_start = bin_info["START"].astype(np.int64).values
    bin_end = bin_info["END"].astype(np.int64).values

    n_bins = len(bin_info)
    n_cl = len(cell_line_cols)
    rt_sums = np.zeros((n_bins, n_cl), dtype=np.float64)
    rt_counts = np.zeros(n_bins, dtype=np.int64)

    chr_to_bin_idx = {}
    for i, c in enumerate(bin_chr):
        chr_to_bin_idx.setdefault(c, []).append(i)

    chr_to_bin_starts = {}
    chr_to_bin_ends = {}
    chr_to_bin_indices = {}
    for c, idx_list in chr_to_bin_idx.items():
        idx_arr = np.array(idx_list, dtype=np.int64)
        chr_to_bin_indices[c] = idx_arr
        chr_to_bin_starts[c] = bin_start[idx_arr]
        chr_to_bin_ends[c] = bin_end[idx_arr]

    logging.info(f"built bin index for {len(chr_to_bin_starts)} chromosomes, {n_bins} bins")

    rt_values = rt_raw[cell_line_cols].values.astype(np.float64)
    total_mapped = 0
    all_mapped_global = []
    all_rt_vals = []
    for chrom, chrom_df in rt_raw.groupby("Chr", sort=False):
        if chrom not in chr_to_bin_starts:
            logging.info(f"{chrom}: skipped (no bins)")
            continue
        starts = chr_to_bin_starts[chrom]
        ends = chr_to_bin_ends[chrom]
        global_idx = chr_to_bin_indices[chrom]

        positions = chrom_df["Position"].values
        bin_i = np.searchsorted(starts, positions, side="right") - 1
        valid = (bin_i >= 0) & (bin_i < len(starts)) & (positions < ends[np.clip(bin_i, 0, len(starts) - 1)])
        n_mapped = int(valid.sum())
        total_mapped += n_mapped
        logging.info(f"{chrom}: {len(chrom_df)} SNPs, {n_mapped} mapped to {len(starts)} bins")
        bin_i = bin_i[valid]
        all_mapped_global.append(global_idx[bin_i])
        all_rt_vals.append(rt_values[chrom_df.index[valid]])

    logging.info(f"total RT SNPs mapped to bins: {total_mapped}/{len(rt_raw)}")

    # Vectorized aggregation with np.bincount (much faster than np.add.at)
    if all_mapped_global:
        all_mapped_global = np.concatenate(all_mapped_global)
        all_rt_vals = np.concatenate(all_rt_vals)
        rt_counts = np.bincount(all_mapped_global, minlength=n_bins).astype(np.int64)
        for ci in range(n_cl):
            rt_sums[:, ci] = np.bincount(
                all_mapped_global, weights=all_rt_vals[:, ci], minlength=n_bins
            )

    with np.errstate(invalid="ignore"):
        rt_means = rt_sums / rt_counts[:, np.newaxis]
    rt_means[rt_counts == 0] = np.nan

    result = bin_info[["#CHR", "START", "END"]].copy()
    for ci, col in enumerate(cell_line_cols):
        result[col] = rt_means[:, ci]

    n_missing = int(np.sum(rt_counts == 0))
    n_total = n_bins
    frac_missing = n_missing / n_total
    if frac_missing > 0.05:
        logging.warning(
            f"RT: {n_missing}/{n_total} bins ({frac_missing:.1%}) have no overlapping "
            f"RT SNPs; these will be filled with NaN and excluded from fitting."
        )
    logging.info(f"RT loaded for {n_cl} cell lines, {n_total} bins")
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
    gc_quantile=(0.05, 0.95),
    eps=1e-8,
):
    """Correct RDR for GC (+ mappability, + RT) bias via natural cubic splines in log2 space.

    Fits OLS on a natural spline basis (df=5) of GC per sample, with optional
    mappability and replication-timing covariates. Corrected RDR is computed
    from residuals shifted to the sample mean, then exponentiated and
    mean-normalized to 1.

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
    eps : float
        Small constant added before log2 to avoid log(0).

    Returns
    -------
    np.ndarray
        Corrected RDR matrix (bins x samples).
    """
    logging.info("ASCAT-style spline RDR correction")
    gc = gc_df["GC"].to_numpy()
    mapv = gc_df["MAP"].to_numpy()

    gc_lo, gc_hi = np.nanquantile(gc, gc_quantile)
    logging.info(f"GC-content range for fitting: [{gc_lo}, {gc_hi}]")

    autosomes_mask = _make_autosomes_mask(gc_df["#CHR"].astype(str).to_numpy())

    pdf = PdfPages(os.path.join(out_dir, "correction_diagnostics.pdf")) if out_dir else None

    corrected_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        raw_rdrs = raw_rdr_mat[:, si]
        log_rdr = np.log2(raw_rdrs + eps)
        rt_vals, best_rt_name = _select_rt_for_sample(raw_rdrs, rt_df, autosomes_mask)

        cov_df = pd.DataFrame({"log_rdr": log_rdr, "GC": gc})
        formula_parts = ["cr(GC, df=5)"]
        if has_mapp:
            cov_df["MAP"] = mapv
            formula_parts.append("cr(MAP, df=5)")
        if rt_vals is not None:
            cov_df["RT"] = rt_vals
            formula_parts.append("cr(RT, df=5)")

        formula_rhs = " + ".join(formula_parts)

        fit_mask = (
            (gc >= gc_lo) & (gc <= gc_hi) &
            np.isfinite(log_rdr)
        )
        if has_mapp:
            fit_mask &= np.isfinite(mapv)
        if rt_vals is not None:
            fit_mask &= np.isfinite(rt_vals)

        logging.info(f"{rep_id}: fitting spline on {fit_mask.sum()}/{len(fit_mask)} bins")

        fit_df = cov_df.loc[fit_mask].copy()
        y_fit = fit_df["log_rdr"].values
        X_fit = patsy.dmatrix(f"{formula_rhs} - 1", data=fit_df, return_type="dataframe")

        ols_model = sm_ols.OLS(y_fit, X_fit).fit()
        logging.info(f"{rep_id}: spline OLS R²={ols_model.rsquared:.4f}")

        pred_df = cov_df.copy()
        for col in pred_df.columns:
            if col == "log_rdr":
                continue
            pred_df[col] = pred_df[col].fillna(pred_df[col].median())
        X_all = patsy.dmatrix(
            f"{formula_rhs} - 1", data=pred_df, return_type="dataframe"
        )
        predicted = ols_model.predict(X_all)

        residuals = log_rdr - predicted
        mean_log_rdr_fit = np.mean(log_rdr[fit_mask])
        corrected_log_rdr = residuals + mean_log_rdr_fit

        corrected_rdr = np.exp2(corrected_log_rdr)

        nan_mask = ~np.isfinite(log_rdr)
        if rt_vals is not None:
            nan_mask |= ~np.isfinite(rt_vals)
        if nan_mask.any():
            med_corr = np.nanmedian(corrected_rdr[~nan_mask])
            corrected_rdr[nan_mask] = med_corr
            logging.info(
                f"{rep_id}: imputed {nan_mask.sum()} NaN bins with median={med_corr:.4f}"
            )

        corrected_mat[:, si] = corrected_rdr / np.mean(corrected_rdr)

        if pdf is not None:
            fitted_rdr = np.exp2(predicted)
            plot_correction_diagnostics(
                gc, raw_rdrs, corrected_mat[:, si],
                fitted_rdr=fitted_rdr, rep_id=rep_id, pdf=pdf,
                rt_vals=rt_vals, rt_name=best_rt_name,
            )

    if pdf is not None:
        pdf.close()
    return corrected_mat
