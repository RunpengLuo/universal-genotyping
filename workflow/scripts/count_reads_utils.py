"""Utility functions for window/bin depth processing and QC plotting.

Used by rd_correct.py, combine_counts.py, and compute_rdr_bulk_bin.py.
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


def log_mad_and_plot(
    pos_df,
    mat,
    labels,
    genome_size,
    qc_dir,
    prefix,
    unit,
    val_type,
    ylim,
    **plot_kwargs,
):
    """Log MAD per column and generate 1-D chromosome scatter plots."""
    for i, label in enumerate(labels):
        v = mat[:, i]
        m = np.isfinite(v)
        if m.any():
            mad = np.median(np.abs(v[m] - np.median(v[m])))
            logging.info(f"  {prefix:<28s} {label:<8s}: MAD={mad:.4f}")
        plot_file = os.path.join(qc_dir, f"{prefix}.{label}.pdf")
        plot_1d_sample(
            pos_df,
            v,
            genome_size,
            plot_file,
            unit=unit,
            val_type=val_type,
            max_ylim=ylim,
            **plot_kwargs,
        )


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
    smooth: bool = False,
    smooth_frac: float = 0.05,
    smooth_color: str = "orange",
    smooth_linewidth: float = 1.5,
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
        if smooth and m.sum() > 10:
            from statsmodels.nonparametric.smoothers_lowess import lowess as _lowess

            sx = x[m]
            sy = y[m]
            order = np.argsort(sx)
            sx, sy = sx[order], sy[order]
            max_pts = 10000
            if len(sx) > max_pts:
                idx = np.linspace(0, len(sx) - 1, max_pts, dtype=int)
                sx_fit, sy_fit = sx[idx], sy[idx]
            else:
                sx_fit, sy_fit = sx, sy
            fit = _lowess(sy_fit, sx_fit, frac=smooth_frac, return_sorted=True)
            ax.plot(
                fit[:, 0],
                fit[:, 1],
                color=smooth_color,
                linewidth=smooth_linewidth,
                zorder=5,
            )
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
