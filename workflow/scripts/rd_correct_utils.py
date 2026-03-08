"""Utility functions for per-window read depth bias correction.

Contains HMMcopy-style LOWESS correction and GC correction QC plots.
"""

import logging

import numpy as np
from scipy.interpolate import interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import pearsonr, spearmanr, gaussian_kde

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def correct_readcount(
    reads,
    gc,
    mappability=None,
    repliseq=None,
    samplesize=50000,
    routlier=0.01,
    doutlier=0.001,
    min_mappability=0.9,
    lowess_frac_tight=0.03,
    lowess_frac_smooth=0.3,
    grid_size=1001,
    seed=42,
):
    """HMMcopy-style correctReadcount: GC + optional mappability + optional
    replication timing stages (ACEseq-style sequential LOWESS correction).

    Parameters
    ----------
    reads : np.ndarray
        Raw read counts per window (1-D).
    gc : np.ndarray
        Per-window GC fraction in [0, 1].
    mappability : np.ndarray or None
        Per-window mappability values in [0, 1]. If None, stage 2 is skipped.
    repliseq : np.ndarray or None
        Per-window consensus replication timing score. If None, stage 3 is
        skipped.  Higher values = earlier replication = higher expected
        coverage in cycling cells.
    samplesize : int
        Max number of ideal bins for loess fitting.
    routlier : float
        Upper quantile for read-count outlier removal.
    doutlier : float
        Quantile for GC/mappability domain outlier removal (top/bottom).
    min_mappability : float
        Minimum mappability for ideal bins.
    lowess_frac_tight : float
        Bandwidth fraction for the first (tight) LOWESS pass.
    lowess_frac_smooth : float
        Bandwidth fraction for the second (smoothing) LOWESS pass.
    grid_size : int
        Number of grid points for LOWESS interpolation.
    seed : int
        Base random seed for subsampling ideal bins (incremented per stage).

    Returns
    -------
    np.ndarray
        Corrected read counts (same length as *reads*).  NaN where
        correction is not possible.
    float
        SSE from the GC LOWESS fit (sum of squared residuals between raw
        depth and LOWESS prediction at valid bins).
    """
    reads = reads.astype(np.float64)
    n = len(reads)

    def _dedup_sorted(xy):
        """Remove duplicate x-values from sorted LOWESS output to avoid
        division-by-zero in interp1d."""
        _, idx = np.unique(xy[:, 0], return_index=True)
        return xy[idx]

    def _dedup_input(y, x):
        """Average y-values for duplicate x-values to avoid degenerate
        local regressions (division by zero) inside LOWESS."""
        ux, inv = np.unique(x, return_inverse=True)
        if len(ux) == len(x):
            return y, x
        uy = np.zeros(len(ux), dtype=np.float64)
        np.add.at(uy, inv, y)
        counts = np.bincount(inv).astype(np.float64)
        uy /= counts
        return uy, ux

    def _fit_lowess_interp(y, x, grid):
        """Tight LOWESS -> grid smooth -> final interpolator."""
        y, x = _dedup_input(y, x)
        if len(x) < 2:
            return None
        s1 = lowess(y, x, frac=lowess_frac_tight, return_sorted=True)
        s1 = _dedup_sorted(s1)
        if len(s1) < 2:
            return None
        s1_fn = interp1d(
            s1[:, 0],
            s1[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        s2 = lowess(s1_fn(grid), grid, frac=lowess_frac_smooth, return_sorted=True)
        s2 = _dedup_sorted(s2)
        if len(s2) < 2:
            return None
        return interp1d(
            s2[:, 0],
            s2[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )

    def _apply_stage(
        prev,
        covariate,
        cov_lo,
        cov_hi,
        grid,
        seed,
        stage_name="",
        extra_ideal_mask=None,
    ):
        """Apply one LOWESS correction stage and return (corrected, sse)."""
        valid = (prev > 0) & np.isfinite(covariate)
        val_hi = np.nanquantile(prev[valid], 1.0 - routlier)
        ideal = valid & (prev <= val_hi) & (covariate >= cov_lo) & (covariate <= cov_hi)
        if extra_ideal_mask is not None:
            ideal &= extra_ideal_mask

        ideal_idx = np.where(ideal)[0]
        if ideal_idx.size < 10:
            logging.warning(
                f"    {stage_name}: only {ideal_idx.size} ideal bins; skipping correction"
            )
            return prev.copy(), 0.0
        if ideal_idx.size > samplesize:
            rng = np.random.default_rng(seed)
            ideal_idx = rng.choice(ideal_idx, size=samplesize, replace=False)

        interp_fn = _fit_lowess_interp(prev[ideal_idx], covariate[ideal_idx], grid)
        if interp_fn is None:
            logging.warning(
                f"    {stage_name}: LOWESS fit failed (too few distinct values); skipping correction"
            )
            return prev.copy(), 0.0
        predicted = interp_fn(covariate)

        # SSE between raw values and LOWESS prediction (valid bins only)
        valid_pred = (predicted > 0) & (prev > 0)
        sse = float(np.sum((prev[valid_pred] - predicted[valid_pred]) ** 2))

        with np.errstate(invalid="ignore", divide="ignore"):
            corrected = np.where(
                (predicted > 0) & (prev > 0),
                prev / predicted,
                np.nan,
            )
        valid_corr = (predicted > 0) & (prev > 0)
        if valid_corr.any():
            scale = np.median(prev[valid_corr]) / np.median(corrected[valid_corr])
            corrected[valid_corr] *= scale

        n_valid = int(valid.sum())
        n_ideal = ideal_idx.size
        n_nan = int(np.isnan(corrected).sum())
        ideal_pct = n_ideal / max(n_valid, 1) * 100
        nan_pct = n_nan / max(n, 1) * 100
        logging.info(
            f"    {stage_name:<5s}  {n_ideal:>8d}/{n_valid} ({ideal_pct:5.1f}%) fit,  "
            f"{n_nan:>8d}/{n} ({nan_pct:5.1f}%) NaN"
        )
        return corrected, sse

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)
    grid_01 = np.linspace(0, 1, grid_size)
    extra = (mappability >= min_mappability) if mappability is not None else None
    cor, gc_sse = _apply_stage(
        reads,
        gc,
        gc_lo,
        gc_hi,
        grid_01,
        seed=seed,
        stage_name="GC",
        extra_ideal_mask=extra,
    )

    if mappability is not None:
        map_lo = np.nanquantile(mappability, doutlier)
        map_hi = np.nanquantile(mappability, 1.0 - doutlier)
        cor, _ = _apply_stage(
            cor,
            mappability,
            map_lo,
            map_hi,
            grid_01,
            seed=seed + 1,
            stage_name="MAP",
        )

    if repliseq is not None:
        repli_finite = np.isfinite(repliseq)
        repli_lo = np.nanquantile(repliseq[repli_finite], doutlier)
        repli_hi = np.nanquantile(repliseq[repli_finite], 1.0 - doutlier)
        grid_repli = np.linspace(repli_lo, repli_hi, grid_size)
        cor, _ = _apply_stage(
            cor,
            repliseq,
            repli_lo,
            repli_hi,
            grid_repli,
            seed=seed + 2,
            stage_name="REPLI",
        )

    return cor.astype(np.float32), gc_sse


def correct_readcount_median(
    reads,
    gc,
    mappability=None,
    repliseq=None,
    doutlier=0.001,
    min_mappability=0.9,
    eps_quantile=0.01,
):
    """Quadratic median quantile regression GC correction.

    Fits RD ~ GC + GC**2 via median quantile regression on valid bins,
    then divides raw depth by predicted and rescales to preserve median.
    Mappability is used as a filter (not a covariate).
    Replication timing correction is not yet implemented (TODO).

    Returns
    -------
    np.ndarray
        Corrected read counts. NaN where correction is not possible.
    float
        SSE from the GC fit.
    """
    import pandas as pd
    import statsmodels.formula.api as smf

    reads = reads.astype(np.float64)
    n = len(reads)

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)

    # Filter: positive depth, GC in range, sufficient mappability
    valid = (reads > 0) & np.isfinite(gc) & (gc >= gc_lo) & (gc <= gc_hi)
    if mappability is not None:
        valid &= mappability >= min_mappability

    if repliseq is not None:
        logging.info("    MEDIAN: RT correction not yet implemented; ignoring repliseq")

    # Build fitting DataFrame — GC-only regression
    fit_df = pd.DataFrame({"RD": reads[valid], "GC": gc[valid]})
    pred_df = pd.DataFrame({"GC": gc})
    formula = "RD ~ GC + I(GC**2)"

    n_fit = len(fit_df)
    logging.info(f"    MEDIAN  {n_fit:>8d}/{n} fitting bins, formula: {formula}")

    if n_fit < 10:
        logging.warning("    MEDIAN: too few fitting bins; skipping correction")
        return reads.astype(np.float32), 0.0

    res = smf.quantreg(formula, data=fit_df).fit(q=0.5)
    predicted = res.predict(pred_df).to_numpy()

    # SSE between raw and predicted at valid bins
    all_valid = (reads > 0) & np.isfinite(gc)
    sse = float(np.sum((reads[all_valid] - predicted[all_valid]) ** 2))

    # Correct: divide raw by predicted, clip predicted floor
    eps = (
        np.nanquantile(predicted[predicted > 0], eps_quantile)
        if (predicted > 0).any()
        else 1.0
    )
    den = np.clip(np.nan_to_num(predicted, nan=eps), eps, None)

    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = np.where(reads > 0, reads / den, np.nan)

    # Set bins below mappability threshold to NaN
    if mappability is not None:
        corrected[mappability < min_mappability] = np.nan

    # Rescale to preserve median of valid bins
    valid_corr = np.isfinite(corrected) & (reads > 0)
    if valid_corr.any():
        scale = np.median(reads[valid_corr]) / np.median(corrected[valid_corr])
        corrected[valid_corr] *= scale

    n_nan = int(np.isnan(corrected).sum())
    logging.info(f"    MEDIAN  {n_nan:>8d}/{n} ({n_nan / max(n, 1) * 100:5.1f}%) NaN")

    return corrected.astype(np.float32), sse


def _rolling_median(x, fraction):
    """Rolling median with mirrored edge padding.

    Reimplements the smoothing approach from CNVkit (cnvlib/smoothing.py).
    """
    import pandas as pd

    x = np.asarray(x, dtype=float)
    wing = max(3, int(np.ceil(len(x) * fraction * 0.5)))
    wing = min(wing, len(x) - 1)
    padded = np.concatenate((x[wing - 1 :: -1], x, x[: -wing - 1 : -1]))
    rolled = pd.Series(padded).rolling(2 * wing + 1, 1, center=True).median()
    return np.asarray(rolled[wing:-wing], dtype=float)


def _center_by_window(log2_vals, sort_key, fraction=None, seed=0xA5EED):
    """CNVkit-style bias correction via rolling median on shuffled+sorted data.

    Reimplements cnvlib/fix.py:center_by_window. The key insight is:
    1. Shuffle bins randomly to break spatial clustering of same-GC bins
       (prevents re-centering of real CNV regions)
    2. Sort by the bias covariate (e.g. GC content)
    3. Rolling median to estimate the bias trend
    4. Subtract bias from log2 values

    Reference: https://cnvkit.readthedocs.io/en/stable/bias.html

    Parameters
    ----------
    log2_vals : np.ndarray
        Log2 depth values.
    sort_key : np.ndarray
        Per-bin bias covariate (e.g. GC fraction).
    fraction : float or None
        Rolling window fraction. If None, uses max(0.01, n^-0.5) following
        the CNVkit default.
    seed : int
        Random seed for shuffling (CNVkit uses 0xA5EED).

    Returns
    -------
    np.ndarray
        Corrected log2 values (same length, original order).
    """
    n = len(log2_vals)
    if n < 10:
        return log2_vals.copy()

    if fraction is None:
        fraction = max(0.01, n**-0.5)

    log2 = log2_vals.copy()

    # 1. Shuffle to break spatial clustering
    rng = np.random.default_rng(seed)
    shuffle_order = rng.permutation(n)
    log2_shuf = log2[shuffle_order]
    key_shuf = sort_key[shuffle_order]

    # 2. Sort by bias covariate
    sort_order = np.argsort(key_shuf, kind="mergesort")
    log2_sorted = log2_shuf[sort_order]

    # 3. Rolling median to estimate bias
    biases = _rolling_median(log2_sorted, fraction)

    # 4. Subtract bias
    log2_sorted -= biases

    # Reverse sort and shuffle to restore original order
    unsort = np.empty_like(sort_order)
    unsort[sort_order] = np.arange(n)
    log2_unshuf = log2_sorted[unsort]

    unshuffle = np.empty_like(shuffle_order)
    unshuffle[shuffle_order] = np.arange(n)
    return log2_unshuf[unshuffle]


def correct_readcount_wes(reads, gc, is_target, mappability=None, min_mappability=0.9):
    """WES-specific bias correction using CNVkit-style rolling median.

    Applies GC correction separately to target and antitarget windows using
    the center_by_window approach from CNVkit (cnvlib/fix.py). Works in log2
    space and corrects for GC bias via rolling median on shuffled+sorted data.

    Optionally applies a second pass for mappability correction.

    Parameters
    ----------
    reads : np.ndarray
        Raw read counts per window (1-D).
    gc : np.ndarray
        Per-window GC fraction in [0, 1].
    is_target : np.ndarray
        Boolean mask indicating target (True) vs antitarget (False) windows.
    mappability : np.ndarray or None
        Per-window mappability values in [0, 1]. If provided, a second
        correction pass is applied for mappability bias.
    min_mappability : float
        Minimum mappability for valid bins (bins below are set to NaN).

    Returns
    -------
    np.ndarray
        Corrected read counts (same length as *reads*).
    float
        SSE from the GC correction (sum of squared residuals between raw
        and corrected depth at valid bins, accumulated across target and
        antitarget groups).
    """
    corrected = np.full(len(reads), np.nan, dtype=np.float64)
    gc_sse = 0.0

    for label, mask in [("target", is_target), ("antitarget", ~is_target)]:
        idx = np.where(mask)[0]
        if len(idx) < 10:
            continue

        r = reads[idx].astype(np.float64)
        g = gc[idx]

        # Valid = positive depth, finite GC, sufficient mappability
        valid = (r > 0) & np.isfinite(g)
        if mappability is not None:
            valid &= mappability[idx] >= min_mappability

        if valid.sum() < 10:
            logging.warning(
                f"    {label}: only {valid.sum()} valid bins; skipping correction"
            )
            continue

        # Convert to log2 space (CNVkit convention)
        with np.errstate(divide="ignore", invalid="ignore"):
            log2_r = np.where(valid, np.log2(np.maximum(r, 1e-10)), np.nan)

        # Center before correction
        median_log2 = np.nanmedian(log2_r[valid])
        log2_centered = log2_r - median_log2

        # GC correction on valid bins
        valid_log2 = log2_centered[valid]
        valid_gc = g[valid]
        corrected_log2 = _center_by_window(valid_log2, valid_gc)

        # Optional mappability correction
        if mappability is not None:
            valid_map = mappability[idx][valid]
            corrected_log2 = _center_by_window(corrected_log2, valid_map)

        # Restore median level and convert back to linear
        corrected_log2 += median_log2
        out = np.full(len(idx), np.nan)
        out[valid] = 2.0**corrected_log2
        corrected[idx] = out

        # SSE: difference between raw and corrected depth at valid bins
        gc_sse += float(np.sum((r[valid] - out[valid]) ** 2))

        n_valid = int(valid.sum())
        n_nan = int(np.isnan(out).sum())
        logging.info(
            f"    {label:<12s}  {n_valid:>8d}/{len(idx)} valid,  "
            f"{n_nan:>8d}/{len(idx)} NaN"
        )

    return corrected.astype(np.float32), gc_sse


def _plot_cov_panel(ax, covariate, reads, xlabel, show_ylabel, rep_id,
                    sse=None, is_before=False):
    """KDE density scatter of readcov vs a covariate on a single axes."""
    valid = (reads > 0) & np.isfinite(covariate)
    if valid.sum() < 20:
        ax.set_visible(False)
        return
    x, y = covariate[valid], reads[valid]
    ylim = np.nanquantile(y, 0.99) * 1.1

    n_pts = len(x)
    rng = np.random.default_rng(0)
    if n_pts > 20000:
        idx = rng.choice(n_pts, size=20000, replace=False)
        kde = gaussian_kde(np.vstack([x[idx], y[idx]]))
    else:
        kde = gaussian_kde(np.vstack([x, y]))

    xgrid = np.linspace(x.min(), x.max(), 200)
    ygrid = np.linspace(0, ylim * 1.5, 200)
    xx, yy = np.meshgrid(xgrid, ygrid)
    positions = np.vstack([xx.ravel(), yy.ravel()])
    zz = kde(positions).reshape(xx.shape)

    ax.pcolormesh(xx, yy, zz, shading="gouraud", cmap="Blues", rasterized=True)
    ax.contour(
        xx, yy, zz, levels=6, colors="steelblue", linewidths=0.5, alpha=0.5
    )

    mad = np.median(np.abs(y - np.median(y)))
    r_pearson, _ = pearsonr(x, y)
    r_spearman, _ = spearmanr(x, y)
    metrics = f"MAD={mad:.2f}  r={r_pearson:.4f}  rho={r_spearman:.4f}"
    if sse is not None and is_before:
        metrics = f"SSE={sse:.1f}  " + metrics
    logging.info(f"  {xlabel:<8s} {rep_id:<8s}: {metrics}")
    ax.set_xlabel(xlabel)
    if show_ylabel:
        ax.set_ylabel("Observed Readcov")
    ax.set_title(f"{rep_id}\n{metrics}", fontsize=8)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, ylim * 1.1)


def plot_gc_correction_pdf(gc, dp_before, dp_after, rep_ids, pdf, gc_sse=None,
                           mappability=None, repliseq=None):
    """Two-page PDF: before/after correction KDE density plots.

    Each page has up to 3 rows (GC, MAP, RT) x nsamples columns.

    Parameters
    ----------
    gc_sse : list of float or None
        Per-sample SSE from the GC fit. Shown on the "Before" panel only.
    mappability : np.ndarray or None
        Per-window mappability values. If provided, a MAP row is added.
    repliseq : np.ndarray or None
        Per-window replication timing values. If provided, an RT row is added.
    """
    nsamples = len(rep_ids)
    panel_w = max(5, 5 * nsamples)

    # Build list of (row_label, covariate, xlabel)
    rows = [("GC", gc, "GC Content")]
    if mappability is not None:
        rows.append(("MAP", mappability, "Mappability"))
    if repliseq is not None:
        rows.append(("RT", repliseq, "Replication Timing"))
    nrows = len(rows)

    is_before = True
    for title, dp_mat in [
        ("Before Correction", dp_before),
        ("After Correction", dp_after),
    ]:
        fig, axes = plt.subplots(
            nrows, nsamples, figsize=(panel_w, 5 * nrows), squeeze=False,
        )
        for ri, (row_label, covariate, xlabel) in enumerate(rows):
            for si, rep_id in enumerate(rep_ids):
                sse = gc_sse[si] if (gc_sse is not None and row_label == "GC") else None
                _plot_cov_panel(
                    axes[ri, si],
                    covariate,
                    dp_mat[:, si],
                    xlabel,
                    show_ylabel=(si == 0),
                    rep_id=rep_id,
                    sse=sse,
                    is_before=is_before,
                )
        fig.suptitle(title, fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        is_before = False
