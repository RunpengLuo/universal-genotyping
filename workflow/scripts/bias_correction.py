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

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def bias_correction_rdr_quantreg(
    raw_rdr_mat: np.ndarray,
    gc_df: pd.DataFrame,
    rep_ids: list,
    has_mapp=False,
    rt_vals=None,
    out_dir=None,
    eps_quantile=0.01,
    gc_quantile=[0.05, 0.95],
):
    """Correct RDR for GC (+ mappability, + RT) bias via median quantile regression.

    Fits a quadratic median quantile regression of RDR on GC per sample,
    with optional mappability and replication-timing covariates. Corrected
    RDR is raw RDR divided by the fitted expected RDR, then mean-normalized
    to 1.

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
    rt_vals : np.ndarray or None
        Per-bin replication timing values (1-D), or None to skip RT.
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

    pdf = PdfPages(os.path.join(out_dir, "correction_diagnostics.pdf")) if out_dir else None

    corrected_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        raw_rdrs = raw_rdr_mat[:, si]

        # Build formula, fit_df, pred_df dynamically
        sample_fit_mask = fit_mask.copy()
        if rt_vals is not None:
            sample_fit_mask &= np.isfinite(rt_vals)

        fit_df = pd.DataFrame({"RD": raw_rdrs[sample_fit_mask], "GC": gc[sample_fit_mask]})
        pred_df = pd.DataFrame({"GC": gc})
        formula = "RD ~ GC + I(GC**2)"
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
        exp_rdrs = res.predict(pred_df).to_numpy()

        eps = np.nanquantile(exp_rdrs, eps_quantile)
        logging.info(f"inferred epsilon for exp RDR={eps}")
        den = np.clip(np.nan_to_num(exp_rdrs, nan=eps), eps, None)
        corr_rdrs = raw_rdrs / den

        corr_factor = np.mean(corr_rdrs)
        logging.info(f"normalization factor={corr_factor}")

        corrected_mat[:, si] = corr_rdrs / corr_factor
        if pdf is not None:
            plot_correction_diagnostics(
                gc, raw_rdrs, corrected_mat[:, si],
                fitted_rdr=exp_rdrs, rep_id=rep_id, pdf=pdf,
                rt_vals=rt_vals, rt_name="REPLI",
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


def gc_correct_depth_loess(
    dp_mat: np.ndarray,
    gc: np.ndarray,
    samplesize: int = 50000,
    routlier: float = 0.01,
    doutlier: float = 0.001,
    mappability: np.ndarray = None,
    min_mappability: float = 0.9,
) -> np.ndarray:
    """Correct raw read depth for GC bias using two-stage loess (Benjamini & Speed 2012).

    Applies per-sample GC correction on raw depth *before* computing RDR.
    The two-stage approach follows HMMcopy: a tight loess on ideal bins,
    then a smoothing pass on a fine GC grid.

    Parameters
    ----------
    dp_mat : np.ndarray
        Raw depth matrix (bins x samples).
    gc : np.ndarray
        Per-bin GC fraction in [0, 1].
    samplesize : int
        Maximum number of ideal bins to use for loess fitting.
    routlier : float
        Upper quantile for read count outlier removal (e.g. 0.01 = top 1%).
    doutlier : float
        Quantile for GC domain outlier removal (top/bottom, e.g. 0.001 = 0.1%).
    mappability : np.ndarray or None
        Per-bin mappability values. If provided, only bins with
        mappability >= *min_mappability* are used for fitting.
    min_mappability : float
        Minimum mappability threshold for ideal bins.

    Returns
    -------
    np.ndarray
        GC-corrected depth matrix (same shape as *dp_mat*).
    """
    n_bins, n_samples = dp_mat.shape
    corrected = np.copy(dp_mat).astype(np.float64)

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)

    for si in range(n_samples):
        reads = dp_mat[:, si].astype(np.float64)

        # Step 1: valid bins (positive reads, valid GC)
        valid = (reads > 0) & np.isfinite(gc) & (gc >= 0)

        # Step 2: ideal bins for fitting
        read_hi = np.nanquantile(reads[valid], 1.0 - routlier)
        ideal = valid & (reads <= read_hi) & (gc >= gc_lo) & (gc <= gc_hi)
        if mappability is not None:
            ideal &= (mappability >= min_mappability)

        ideal_idx = np.where(ideal)[0]
        logging.info(
            f"sample {si}: {ideal_idx.size} ideal bins out of {int(valid.sum())} valid"
        )

        # Step 3: subsample
        if ideal_idx.size > samplesize:
            rng = np.random.default_rng(42 + si)
            ideal_idx = rng.choice(ideal_idx, size=samplesize, replace=False)

        gc_ideal = gc[ideal_idx]
        reads_ideal = reads[ideal_idx]

        # Step 4a: stage 1 — tight loess on ideal bins
        stage1 = lowess(reads_ideal, gc_ideal, frac=0.03, return_sorted=True)
        # stage1 is (N, 2) sorted by gc; columns are (gc, fitted)

        # Step 4b: stage 2 — evaluate stage-1 on fine grid, then smooth
        grid = np.linspace(0, 1, 1001)
        stage1_interp = interp1d(
            stage1[:, 0], stage1[:, 1],
            kind="linear", bounds_error=False, fill_value="extrapolate",
        )
        grid_pred = stage1_interp(grid)
        stage2 = lowess(grid_pred, grid, frac=0.3, return_sorted=True)
        final_interp = interp1d(
            stage2[:, 0], stage2[:, 1],
            kind="linear", bounds_error=False, fill_value="extrapolate",
        )

        # Step 5: correct all bins
        predicted = final_interp(gc)
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = np.where(
                (predicted > 0) & (reads > 0),
                reads / predicted,
                0.0,
            )
        # Rescale so corrected median matches original median (for valid bins)
        valid_corr = (predicted > 0) & (reads > 0)
        if valid_corr.any():
            scale = np.median(reads[valid_corr]) / np.median(corr[valid_corr])
            corr[valid_corr] *= scale

        corrected[:, si] = corr

    return corrected.astype(np.float32)


def plot_depth_gc_correction(
    gc, dp_before, dp_after, rep_id, pdf,
):
    """One-page diagnostic: depth vs GC before and after loess correction."""
    mad = lambda x: np.median(np.abs(x - np.median(x)))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    valid = (dp_before > 0) & np.isfinite(gc)
    mad_before = mad(dp_before[valid])
    r_before, _ = pearsonr(gc[valid], dp_before[valid])
    axes[0].scatter(gc[valid], dp_before[valid], s=2, alpha=0.2, rasterized=True)
    axes[0].set_xlabel("GC")
    axes[0].set_ylabel("Read depth")
    axes[0].set_title(f"Before GC correction\nMAD={mad_before:.4f}  r(GC,depth)={r_before:.4f}")

    valid2 = (dp_after > 0) & np.isfinite(gc)
    mad_after = mad(dp_after[valid2])
    r_after, _ = pearsonr(gc[valid2], dp_after[valid2])
    axes[1].scatter(gc[valid2], dp_after[valid2], s=2, alpha=0.2, rasterized=True)
    axes[1].set_xlabel("GC")
    axes[1].set_ylabel("Read depth")
    axes[1].set_title(f"After GC correction\nMAD={mad_after:.4f}  r(GC,depth)={r_after:.4f}")

    fig.suptitle(f"{rep_id} — loess depth correction", fontsize=14)
    plt.tight_layout()
    pdf.savefig(fig, dpi=150)
    plt.close(fig)
    logging.info(
        f"{rep_id}: MAD {mad_before:.4f} -> {mad_after:.4f}, "
        f"r(GC,depth) {r_before:.4f} -> {r_after:.4f}"
    )


def bias_correction_rdr_spline(
    raw_rdr_mat: np.ndarray,
    gc_df: pd.DataFrame,
    rep_ids: list,
    has_mapp=False,
    rt_vals=None,
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
    rt_vals : np.ndarray or None
        Per-bin replication timing values (1-D), or None to skip RT.
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

    pdf = PdfPages(os.path.join(out_dir, "correction_diagnostics.pdf")) if out_dir else None

    corrected_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        raw_rdrs = raw_rdr_mat[:, si]
        log_rdr = np.log2(raw_rdrs + eps)

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
        corr_log_rdrs = residuals + mean_log_rdr_fit

        corr_rdrs = np.exp2(corr_log_rdrs)

        nan_mask = ~np.isfinite(log_rdr)
        if rt_vals is not None:
            nan_mask |= ~np.isfinite(rt_vals)
        if nan_mask.any():
            med_corr = np.nanmedian(corr_rdrs[~nan_mask])
            corr_rdrs[nan_mask] = med_corr
            logging.info(
                f"{rep_id}: imputed {nan_mask.sum()} NaN bins with median={med_corr:.4f}"
            )

        corr_factor = np.mean(corr_rdrs)
        logging.info(f"normalization factor={corr_factor}")

        corrected_mat[:, si] = corr_rdrs / corr_factor

        if pdf is not None:
            fitted_rdr = np.exp2(predicted)
            plot_correction_diagnostics(
                gc, raw_rdrs, corrected_mat[:, si],
                fitted_rdr=fitted_rdr, rep_id=rep_id, pdf=pdf,
                rt_vals=rt_vals, rt_name="REPLI",
            )

    if pdf is not None:
        pdf.close()
    return corrected_mat
