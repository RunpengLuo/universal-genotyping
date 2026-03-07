"""Utility functions for window/bin depth processing and QC plotting.

Used by rd_correct.py and combine_counts.py.
"""

import os
import logging

import numpy as np
import pandas as pd

import pyranges as pr

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import get_chr_sizes, read_region_file


def compute_overlap_weights(win_s, win_e, bb_starts, bb_ends):
    """Find all (window, bin) overlaps and return bp-length weights.

    Returns
    -------
    ov_win : np.ndarray[intp]   window-local indices
    ov_bin : np.ndarray[intp]   bin-local indices
    ov_wt  : np.ndarray[float64] overlap length in bp
    """
    n_bb = len(bb_starts)
    first_bin = np.searchsorted(bb_ends, win_s, side="right")
    last_bin = np.searchsorted(bb_starts, win_e, side="left") - 1

    ov_win_list, ov_bin_list, ov_wt_list = [], [], []
    for wi in range(len(win_s)):
        fb, lb = int(first_bin[wi]), int(last_bin[wi])
        if fb > lb or fb >= n_bb or lb < 0:
            continue
        fb = max(fb, 0)
        lb = min(lb, n_bb - 1)
        ws, we = int(win_s[wi]), int(win_e[wi])
        if we <= ws:
            continue
        for bi in range(fb, lb + 1):
            ov_start = max(ws, int(bb_starts[bi]))
            ov_end = min(we, int(bb_ends[bi]))
            ov_len = ov_end - ov_start
            if ov_len > 0:
                ov_win_list.append(wi)
                ov_bin_list.append(bi)
                ov_wt_list.append(ov_len)

    if not ov_win_list:
        return (
            np.array([], dtype=np.intp),
            np.array([], dtype=np.intp),
            np.array([], dtype=np.float64),
        )
    return (
        np.array(ov_win_list, dtype=np.intp),
        np.array(ov_bin_list, dtype=np.intp),
        np.array(ov_wt_list, dtype=np.float64),
    )


def weighted_bincount_mean(values, bin_idx, weights, n_bins):
    """Weighted mean per bin via bincount.  NaN-aware: ignores NaN values."""
    finite = np.isfinite(values)
    bi = bin_idx[finite]
    wt = weights[finite]
    v = values[finite]
    sums = np.bincount(bi, weights=v * wt, minlength=n_bins)
    wt_sums = np.bincount(bi, weights=wt, minlength=n_bins)
    with np.errstate(invalid="ignore"):
        return np.where(wt_sums > 0, sums / wt_sums, np.nan)


def get_mask_by_region_intervals(
    bins_df: pd.DataFrame, regions: pd.DataFrame
) -> np.ndarray:
    """Return a boolean mask indicating whether each bin overlaps any region.

    Parameters
    ----------
    bins_df : pd.DataFrame
        Must have ``#CHR``, ``START``, ``END`` columns (0-based half-open).
    regions : pd.DataFrame
        Must have ``Chromosome``, ``Start``, ``End`` columns (from ``read_region_file``).
    """
    tmp = bins_df[["#CHR", "START", "END"]].copy()
    tmp["_idx"] = np.arange(len(tmp))
    pr_bins = pr.PyRanges(
        tmp.rename(columns={"#CHR": "Chromosome", "START": "Start", "END": "End"})
    )
    pr_regions = pr.PyRanges(regions)
    overlapping = pr_bins.overlap(pr_regions).df
    keep_idx = set(overlapping["_idx"].unique()) if not overlapping.empty else set()
    return np.isin(np.arange(len(bins_df)), list(keep_idx))


def log_nan_summary(name, mat, labels, n_total):
    """Log per-column NaN counts for a matrix."""
    for i, label in enumerate(labels):
        col = mat[:, i] if mat.ndim == 2 else mat
        n_nan = int(np.isnan(col).sum())
        pct = n_nan / max(n_total, 1) * 100
        logging.info(
            f"  {name:<16s} {label:<8s}: {n_nan:>8d}/{n_total} ({pct:5.1f}%) NaN"
        )


def compute_gc_rd_stats(mat, gc_vals, labels, n_gc_bins=100):
    """Compute per-label Pearson/Spearman corr(RD, GC) and std of binned median RD.

    Parameters
    ----------
    mat : np.ndarray
        (n_bins, n_samples) depth matrix.
    gc_vals : np.ndarray
        Per-bin GC fraction (same length as mat rows).
    labels : list[str]
        Column labels (sample/rep IDs).
    n_gc_bins : int
        Number of equal-width GC bins in [0, 1].

    Returns
    -------
    gc_corr : dict[str, tuple[float, float]]
        {label: (pearson_r, spearman_r)}
    gc_bin_median_std : dict[str, float]
        {label: std of per-GC-bin median RD (A_GC)}
    """
    from scipy.stats import pearsonr, spearmanr

    gc_bins = np.linspace(0, 1, n_gc_bins + 1)
    gc_corr = {}
    gc_bin_median_std = {}

    for i, label in enumerate(labels):
        v = mat[:, i] if mat.ndim == 2 else mat
        valid = np.isfinite(v) & np.isfinite(gc_vals)

        if valid.sum() > 2:
            pr_val, _ = pearsonr(gc_vals[valid], v[valid])
            sr_val, _ = spearmanr(gc_vals[valid], v[valid])
            gc_corr[label] = (pr_val, sr_val)
        else:
            gc_corr[label] = (np.nan, np.nan)

        bin_idx = np.digitize(gc_vals, gc_bins) - 1
        bin_idx = np.clip(bin_idx, 0, n_gc_bins - 1)
        medians = []
        for b in range(n_gc_bins):
            mask = (bin_idx == b) & np.isfinite(v)
            if mask.any():
                medians.append(np.median(v[mask]))
        gc_bin_median_std[label] = float(np.std(medians)) if medians else np.nan

    return gc_corr, gc_bin_median_std


def plot_rd_gc(
    pos_df,
    mat,
    labels,
    genome_size,
    qc_dir,
    prefix,
    unit,
    val_type,
    ylim,
    gc_corr=None,
    gc_bin_median_std=None,
    **plot_kwargs,
):
    """Log GC-correlation stats per column and generate a combined multi-sample
    1-D chromosome scatter PDF.

    Produces a single ``{prefix}.pdf`` with one page per chromosome.  Each page
    contains vertically stacked subplots (one row per sample, shared x-axis).

    Parameters
    ----------
    gc_corr : dict[str, tuple[float, float]] | None
        {label: (pearson_r, spearman_r)} from ``compute_gc_rd_stats``.
    gc_bin_median_std : dict[str, float] | None
        {label: std of per-GC-bin median RD} from ``compute_gc_rd_stats``.
    """
    for i, label in enumerate(labels):
        if gc_corr is not None and label in gc_corr:
            pr_val, sr_val = gc_corr[label]
            logging.info(
                f"  {prefix:<28s} {label:<8s}: "
                f"Pearson={pr_val:.4f}, Spearman={sr_val:.4f}"
            )
        if gc_bin_median_std is not None and label in gc_bin_median_std:
            logging.info(
                f"  {prefix:<28s} {label:<8s}: A_GC std={gc_bin_median_std[label]:.4f}"
            )
    plot_file = os.path.join(qc_dir, f"{prefix}.pdf")
    plot_1d_multi_sample(
        pos_df,
        mat,
        labels,
        genome_size,
        plot_file,
        unit=unit,
        val_type=val_type,
        max_ylim=ylim,
        **plot_kwargs,
    )


def plot_1d_multi_sample(
    pos_df: pd.DataFrame,
    mat: np.ndarray,
    labels: list,
    genome_size: str,
    out_file: str,
    unit="window",
    val_type="RD",
    s=4,
    dpi=72,
    alpha=0.6,
    min_ylim=0.0,
    max_ylim=None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
):
    """Multi-sample 1-D chromosome scatter plot: one page per chromosome,
    vertically stacked subplots (one row per sample, shared x-axis).

    Parameters
    ----------
    mat : np.ndarray
        (n_windows, n_samples) value matrix.
    labels : list[str]
        Sample labels, length == mat.shape[1].
    """
    n_samples = len(labels)
    logging.info(
        f"chromosome wide {unit}-level {val_type} multi-sample plot "
        f"({n_samples} samples), out_file={out_file}"
    )
    chrom_sizes = get_chr_sizes(genome_size)

    _region_by_chr = {}
    if region_bed is not None:
        _region_df = read_region_file(region_bed)
        for _ch, _grp in _region_df.groupby("#CHR", sort=False):
            _region_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    _blacklist_by_chr = {}
    if blacklist_bed is not None:
        _bl_df = read_region_file(blacklist_bed)
        for _ch, _grp in _bl_df.groupby("#CHR", sort=False):
            _blacklist_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    ch = pos_df["#CHR"].to_numpy()
    if "POS" in pos_df.columns:
        pos = pos_df["POS"].to_numpy()
    else:
        pos = ((pos_df["START"].to_numpy() + pos_df["END"].to_numpy()) // 2).astype(
            np.int64
        )
    num_bins = len(pos_df)
    change = np.flatnonzero(ch[1:] != ch[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, num_bins]
    chroms = ch[starts]

    pdf_fd = PdfPages(out_file)
    for chrom, lo, hi in zip(chroms, starts, ends):
        chr_end = chrom_sizes.get(chrom)
        if chr_end is None:
            logging.warning(f"{chrom}: not found in {genome_size}")
            continue

        x = pos[lo:hi]

        fig, axes = plt.subplots(
            nrows=n_samples,
            ncols=1,
            figsize=(40, 3 * n_samples),
            sharex=True,
            squeeze=False,
        )
        axes = axes[:, 0]

        for si, (ax, label) in enumerate(zip(axes, labels)):
            y = mat[lo:hi, si] if mat.ndim == 2 else mat[lo:hi]
            m = np.isfinite(y)

            if chrom in _region_by_chr:
                for reg_start, reg_end in _region_by_chr[chrom]:
                    ax.axvspan(
                        reg_start, reg_end, color="lightblue", alpha=0.15, zorder=0
                    )
            if chrom in _blacklist_by_chr:
                for bl_start, bl_end in _blacklist_by_chr[chrom]:
                    ax.axvspan(
                        bl_start, bl_end, color="lightcoral", alpha=0.2, zorder=0
                    )

            if m.any():
                ax.scatter(x[m], y[m], s=s, alpha=alpha, rasterized=True)

            if val_type in ["AF", "BAF"]:
                ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
                ax.set_ylim(-0.05, 1.05)
            elif max_ylim is not None:
                ax.set_ylim(min_ylim, max_ylim)

            ax.set_ylabel(label, fontsize=10)
            ax.set_xlim(0, chr_end)
            ax.grid(alpha=0.2)

            if si == 0:
                ax.set_title(f"{val_type} plot - {chrom}")

        axes[-1].set_xlabel(f"Position ({unit})")
        fig.tight_layout()
        pdf_fd.savefig(fig, dpi=dpi)
        plt.close(fig)
    pdf_fd.close()


def plot_1d_sample(
    pos_df: pd.DataFrame,
    val: np.ndarray,
    genome_size: str,
    out_file: str,
    unit="SNP",
    val_type="BAF",
    s=4,
    dpi=72,
    alpha=0.6,
    figsize=(40, 3),
    min_ylim=0.0,
    max_ylim=None,
    mask: np.ndarray | None = None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
):
    """
    plot any features like AF, RDR, etc., in 1d chromosome scatter plot.
    """
    logging.info(f"chromosome wide {unit}-level {val_type} plot, out_file={out_file}")
    chrom_sizes = get_chr_sizes(genome_size)

    _region_by_chr = {}
    if region_bed is not None:
        _region_df = read_region_file(region_bed)
        for _ch, _grp in _region_df.groupby("#CHR", sort=False):
            _region_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    _blacklist_by_chr = {}
    if blacklist_bed is not None:
        _bl_df = read_region_file(blacklist_bed)
        for _ch, _grp in _bl_df.groupby("#CHR", sort=False):
            _blacklist_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    ch = pos_df["#CHR"].to_numpy()
    if "POS" in pos_df.columns:
        pos = pos_df["POS"].to_numpy()
    else:
        pos = ((pos_df["START"].to_numpy() + pos_df["END"].to_numpy()) // 2).astype(
            np.int64
        )
    num_bins = len(pos_df)
    change = np.flatnonzero(ch[1:] != ch[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, num_bins]
    chroms = ch[starts]

    pdf_fd = PdfPages(out_file)
    for chrom, lo, hi in zip(chroms, starts, ends):
        chr_end = chrom_sizes.get(chrom)
        if chr_end is None:
            logging.warning(f"{chrom}: not found in {genome_size}")
            continue

        x = pos[lo:hi]
        y = val[lo:hi]
        m = np.isfinite(y)
        mask_chr = mask[lo:hi] if mask is not None else None

        n_plot = int(m.sum())
        if n_plot == 0:
            logging.warning(f"{chrom}: all {val_type} values are non-finite")
            continue
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        if chrom in _region_by_chr:
            for reg_start, reg_end in _region_by_chr[chrom]:
                ax.axvspan(reg_start, reg_end, color="lightblue", alpha=0.15, zorder=0)
        if chrom in _blacklist_by_chr:
            for bl_start, bl_end in _blacklist_by_chr[chrom]:
                ax.axvspan(bl_start, bl_end, color="lightcoral", alpha=0.2, zorder=0)
        if mask_chr is not None:
            kept = m & mask_chr
            filt = m & ~mask_chr
            if filt.any():
                ax.scatter(
                    x[filt],
                    y[filt],
                    s=s,
                    alpha=0.8,
                    color="red",
                    rasterized=True,
                    label=f"filtered ({filt.sum()})",
                )
            if kept.any():
                ax.scatter(
                    x[kept],
                    y[kept],
                    s=s,
                    alpha=0.8,
                    color="blue",
                    rasterized=True,
                    label=f"kept ({kept.sum()})",
                )
            ax.legend(loc="upper right", fontsize=8, markerscale=2)
        else:
            ax.scatter(x[m], y[m], s=s, alpha=alpha, rasterized=True)
        if val_type in ["AF", "BAF"]:
            ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
            ax.set_ylim(-0.05, 1.05)
        elif max_ylim is not None:
            ax.set_ylim(min_ylim, max_ylim)

        ax.set_xlabel(val_type)
        ax.set_xlim(0, chr_end)
        ax.grid(alpha=0.2)
        ax.set_title(f"{val_type} plot - {chrom}")
        pdf_fd.savefig(fig, dpi=dpi)
        plt.close(fig)
    pdf_fd.close()
    return
