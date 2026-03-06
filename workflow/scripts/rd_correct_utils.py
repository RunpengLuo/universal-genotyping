"""Utility functions for per-window read depth bias correction.

Contains HMMcopy-style LOWESS correction and GC correction QC plots.
"""

import logging

import numpy as np
from scipy.interpolate import interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import pearsonr, gaussian_kde

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
    """
    reads = reads.astype(np.float64)
    n = len(reads)

    def _fit_lowess_interp(y, x, grid):
        """Tight LOWESS -> grid smooth -> final interpolator."""
        s1 = lowess(y, x, frac=lowess_frac_tight, return_sorted=True)
        s1_fn = interp1d(
            s1[:, 0],
            s1[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        s2 = lowess(s1_fn(grid), grid, frac=lowess_frac_smooth, return_sorted=True)
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
        """Apply one LOWESS correction stage and return corrected values."""
        valid = (prev > 0) & np.isfinite(covariate)
        val_hi = np.nanquantile(prev[valid], 1.0 - routlier)
        ideal = valid & (prev <= val_hi) & (covariate >= cov_lo) & (covariate <= cov_hi)
        if extra_ideal_mask is not None:
            ideal &= extra_ideal_mask

        ideal_idx = np.where(ideal)[0]
        if ideal_idx.size > samplesize:
            rng = np.random.default_rng(seed)
            ideal_idx = rng.choice(ideal_idx, size=samplesize, replace=False)

        interp_fn = _fit_lowess_interp(prev[ideal_idx], covariate[ideal_idx], grid)
        predicted = interp_fn(covariate)

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
        return corrected

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)
    grid_01 = np.linspace(0, 1, grid_size)
    extra = (mappability >= min_mappability) if mappability is not None else None
    cor = _apply_stage(
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
        cor = _apply_stage(
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
        cor = _apply_stage(
            cor,
            repliseq,
            repli_lo,
            repli_hi,
            grid_repli,
            seed=seed + 2,
            stage_name="REPLI",
        )

    return cor.astype(np.float32)


def correct_readcount_wes(reads, gc, is_target, **kwargs):
    """WES-specific bias correction (e.g. cnvkit-style median-centering).

    TODO: Implement a WES-appropriate correction method. LOWESS is unsuitable
    for WES because target and antitarget windows have very different size
    distributions and depth profiles. Consider cnvkit-style approaches:
    - Separate median-centering for target vs antitarget
    - Reference-based correction using a panel of normals
    - Edge-effect and length-bias correction for exon captures

    Parameters
    ----------
    reads : np.ndarray
        Raw read counts per window (1-D).
    gc : np.ndarray
        Per-window GC fraction in [0, 1].
    is_target : np.ndarray
        Boolean mask indicating target (True) vs antitarget (False) windows.
    **kwargs
        Additional parameters (TBD).

    Returns
    -------
    np.ndarray
        Corrected read counts (same length as *reads*).
    """
    raise NotImplementedError("correct_readcount_wes is not yet implemented")


def plot_gc_correction_pdf(gc, dp_before, dp_after, rep_ids, pdf):
    """Two-page PDF: before/after GC correction KDE density plots."""
    nsamples = len(rep_ids)
    panel_w = max(5, 5 * nsamples)

    for title, dp_mat in [
        ("Before GC Correction", dp_before),
        ("After GC Correction", dp_after),
    ]:
        fig, axes = plt.subplots(1, nsamples, figsize=(panel_w, 5), squeeze=False)
        for si, rep_id in enumerate(rep_ids):
            ax = axes[0, si]
            reads = dp_mat[:, si]
            valid = (reads > 0) & np.isfinite(gc)
            x, y = gc[valid], reads[valid]
            ylim = np.nanquantile(y, 0.99)

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
            r, _ = pearsonr(x, y)
            logging.info(f"  {title:<24s} {rep_id:<8s}: MAD={mad:.4f}  r={r:.4f}")
            ax.set_xlabel("GC Content")
            if si == 0:
                ax.set_ylabel("Observed Readcov")
            ax.set_title(f"{rep_id}\nMAD={mad:.4f}  r={r:.4f}")
            ax.set_xlim(x.min(), x.max())
            ax.set_ylim(0, ylim * 1.1)
        fig.suptitle(title, fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
