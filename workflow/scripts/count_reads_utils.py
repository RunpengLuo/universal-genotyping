"""Utility functions for window/bin depth processing.

Used by rd_correct.py and combine_counts.py.
"""

import logging

import numpy as np
import pandas as pd

import pyranges as pr



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


