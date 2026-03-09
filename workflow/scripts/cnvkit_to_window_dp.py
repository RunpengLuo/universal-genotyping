"""Convert CNVkit outputs to window.dp.npz + window.tsv.gz.

Uses reference.cnn as the canonical window set. Large windows are subdivided,
blacklisted regions removed, and target/antitarget annotated. Depth from .cnr
is left-joined onto the canonical windows (unmatched → 0.0). Raw .cnn coverage
is loaded for before/after QC comparison when available.
"""

import os
import logging
from snakemake.script import snakemake as sm

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd
import pyranges as pr

from utils import setup_logging, maybe_path, stamp_path
from count_reads_utils import (
    log_nan_summary,
    compute_gc_rd_stats,
    plot_rd_gc,
)
from rd_correct_utils import plot_gc_correction_pdf
from io_utils import read_region_file, get_chr_sizes

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages

setup_logging(sm.log[0])

CHROM_ORDER = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

sample_file = sm.input["sample_file"]
region_bed = sm.input["region_bed"] or None
blacklist_bed = maybe_path(sm.input.get("blacklist_bed", None))
genome_size = sm.input["genome_size"]

cnvkit_dir = sm.params["cnvkit_dir"]
chromosomes = sm.params["chromosomes"]
min_window_length = int(sm.params.get("min_window_length", 1000))

qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
run_id = getattr(sm.params, "run_id", "")

sample_df = pd.read_table(sample_file, sep="\t")
rep_ids = sample_df["REP_ID"].astype(str).tolist()
nsamples = len(sample_df)
target_chroms = {f"chr{c}" for c in chromosomes}

logging.info(f"cnvkit_to_window_dp: {nsamples} samples, {len(target_chroms)} chroms")

ref_cnn = pd.read_table(sm.input["reference_cnn"], sep="\t")
ref_cnn = ref_cnn[ref_cnn["chromosome"].isin(target_chroms)].reset_index(drop=True)

coords = ref_cnn[["chromosome", "start", "end"]].apply(tuple, axis=1)
assert coords.is_unique, "Duplicate coordinates in reference.cnn"

win_df = pd.DataFrame({
    "#CHR": ref_cnn["chromosome"].values,
    "START": ref_cnn["start"].values,
    "END": ref_cnn["end"].values,
})

if "gc" in ref_cnn.columns:
    win_df["GC"] = ref_cnn["gc"].values
    logging.info(f"GC loaded from reference.cnn for {len(win_df)} windows")
else:
    win_df["GC"] = np.nan
    logging.warning("no gc column in reference.cnn")

if "gene" in ref_cnn.columns:
    win_df["is_target"] = (
        ~ref_cnn["gene"].isin(["Antitarget", "Background"])
    ).astype(np.int8)
else:
    win_df["is_target"] = np.int8(1)
    logging.warning("no gene column in reference.cnn; assuming all target")

n_ref = len(win_df)
n_target = int(win_df["is_target"].sum())
logging.info(
    f"reference.cnn: {n_ref} windows ({n_target} target, {n_ref - n_target} antitarget)"
)

def _build_depth_map(df):
    """Build (chromosome, start, end) → depth dict."""
    keys = zip(df["chromosome"], df["start"].astype(int), df["end"].astype(int))
    return dict(zip(keys, df["depth"].astype(float)))


def _fill_depth_matrix(depth_maps, win, fill=0.0):
    """Fill (n_windows, n_samples) matrix from depth maps keyed by coordinates."""
    keys = list(zip(win["#CHR"], win["START"].astype(int), win["END"].astype(int)))
    mat = np.full((len(keys), len(depth_maps)), fill, dtype=np.float32)
    for i, dmap in enumerate(depth_maps):
        mat[:, i] = [dmap.get(k, fill) for k in keys]
    return mat


def _extend_windows_to_midpoints(win_df, blacklist_bed, region_bed, genome_size):
    """Extend each window's START/END to the midpoint of the gap with its neighbor.

    Extensions respect blacklist regions and capture-panel (region) boundaries.
    Edge windows are extended toward chromosome boundaries (0 / chrom_size).
    Windows are never shrunk — only expanded into gaps.
    """
    chr_sizes = get_chr_sizes(genome_size)
    win_df = win_df.copy()

    # Build sorted barrier arrays per chromosome.
    # Left barriers (clamp left extension): blacklist ENDs + all region boundaries.
    # Right barriers (clamp right extension): blacklist STARTs + all region boundaries.
    left_bar_by_chr = {}
    right_bar_by_chr = {}

    bl_by_chr = {}
    if blacklist_bed is not None:
        bl_df = pr.read_bed(blacklist_bed, as_df=True)
        for chrom, grp in bl_df.groupby("Chromosome"):
            bl_by_chr[chrom] = grp[["Start", "End"]].values

    reg_bounds_by_chr = {}
    if region_bed is not None:
        reg_df = read_region_file(region_bed)
        for chrom, grp in reg_df.groupby("#CHR"):
            reg_bounds_by_chr[chrom] = np.unique(
                np.concatenate([grp["START"].values, grp["END"].values])
            )

    for chrom in win_df["#CHR"].unique():
        lb, rb = [], []
        if chrom in bl_by_chr:
            bl = bl_by_chr[chrom]
            lb.append(bl[:, 1])   # blacklist ends
            rb.append(bl[:, 0])   # blacklist starts
        if chrom in reg_bounds_by_chr:
            lb.append(reg_bounds_by_chr[chrom])
            rb.append(reg_bounds_by_chr[chrom])
        left_bar_by_chr[chrom] = np.unique(np.concatenate(lb)) if lb else np.array([], dtype=np.int64)
        right_bar_by_chr[chrom] = np.unique(np.concatenate(rb)) if rb else np.array([], dtype=np.int64)

    starts = win_df["START"].values.copy()
    ends = win_df["END"].values.copy()

    for chrom, idx in win_df.groupby("#CHR", sort=False).groups.items():
        idx = idx.values
        order = np.argsort(starts[idx])
        idx = idx[order]
        n = len(idx)
        if n == 0:
            continue

        s = starts[idx]
        e = ends[idx]
        chrom_size = chr_sizes.get(chrom, int(e[-1]))

        # Midpoints between adjacent windows (vectorized)
        mids = (e[:-1] + s[1:]) // 2
        raw_lefts = np.empty(n, dtype=np.int64)
        raw_rights = np.empty(n, dtype=np.int64)
        raw_lefts[0] = 0
        raw_lefts[1:] = mids
        raw_rights[:-1] = mids
        raw_rights[-1] = chrom_size

        # Clamp left extensions: find max barrier in (raw_lefts[k], s[k]]
        lb = left_bar_by_chr.get(chrom, np.array([], dtype=np.int64))
        if len(lb) > 0:
            lo = np.searchsorted(lb, raw_lefts, side="right")   # first index > raw_left
            hi = np.searchsorted(lb, s, side="right")           # first index > s
            has = hi > lo
            best = lb[np.clip(hi - 1, 0, len(lb) - 1)]
            raw_lefts = np.where(has, np.maximum(raw_lefts, best), raw_lefts)

        # Clamp right extensions: find min barrier in [e[k], raw_rights[k])
        rb = right_bar_by_chr.get(chrom, np.array([], dtype=np.int64))
        if len(rb) > 0:
            lo = np.searchsorted(rb, e, side="left")            # first index >= e
            hi = np.searchsorted(rb, raw_rights, side="left")   # first index >= raw_right
            has = hi > lo
            best = rb[np.clip(lo, 0, len(rb) - 1)]
            raw_rights = np.where(has, np.minimum(raw_rights, best), raw_rights)

        # Never shrink
        starts[idx] = np.minimum(raw_lefts, s)
        ends[idx] = np.maximum(raw_rights, e)

    win_df["START"] = starts
    win_df["END"] = ends
    return win_df


win_keys = set(zip(win_df["#CHR"], win_df["START"].astype(int), win_df["END"].astype(int)))

cnr_depth_maps = []
for rep_id in rep_ids:
    cnr_file = os.path.join(cnvkit_dir, f"{rep_id}.cnr")
    df = pd.read_table(cnr_file, sep="\t")
    df = df[df["chromosome"].isin(target_chroms)].reset_index(drop=True)
    dmap = _build_depth_map(df)
    cnr_depth_maps.append(dmap)
    n_matched = len(win_keys & dmap.keys())
    logging.info(
        f"  {rep_id}: {len(df)} .cnr bins, {n_matched}/{n_ref} matched to reference.cnn"
    )

dp_corr_ref = _fill_depth_matrix(cnr_depth_maps, win_df, fill=0.0)

dp_raw_ref = None
try:
    raw_depth_maps = []
    for rep_id in rep_ids:
        tgt_file = os.path.join(cnvkit_dir, f"{rep_id}.targetcoverage.cnn")
        anti_file = os.path.join(cnvkit_dir, f"{rep_id}.antitargetcoverage.cnn")
        tgt_df = pd.read_table(tgt_file, sep="\t")
        anti_df = pd.read_table(anti_file, sep="\t")
        raw_df = pd.concat([tgt_df, anti_df], ignore_index=True)
        raw_df = raw_df[raw_df["chromosome"].isin(target_chroms)].reset_index(drop=True)
        raw_depth_maps.append(_build_depth_map(raw_df))

    dp_raw_ref = _fill_depth_matrix(raw_depth_maps, win_df, fill=np.nan)
    n_matched = int(np.isfinite(dp_raw_ref[:, 0]).sum())
    logging.info(
        f"raw coverage: {n_matched}/{n_ref} "
        f"({n_matched / max(n_ref, 1) * 100:.1f}%) bins matched"
    )
    if n_matched == 0:
        dp_raw_ref = None
except Exception as e:
    logging.warning(f"could not load raw .cnn files: {e}; skipping before/after comparison")
    dp_raw_ref = None

gc_vals = win_df["GC"].to_numpy()
has_gc = np.isfinite(gc_vals).any()
is_target_mask = win_df["is_target"].to_numpy().astype(bool)

if has_gc:
    def _qc_plot_subset(dp_mat, mask, prefix):
        sub_win = win_df.loc[mask].reset_index(drop=True)
        sub_dp = dp_mat[mask]
        sub_gc = gc_vals[mask]
        if len(sub_win) == 0:
            logging.info(f"  {prefix}: no windows, skipping")
            return None, None
        ylim = max(np.nanquantile(sub_dp, 0.99), 1.0) * 1.1
        min_ylim = min(np.nanquantile(sub_dp, 0.01), 0.0) * 1.1
        gc_corr, gc_std = compute_gc_rd_stats(sub_dp, sub_gc, rep_ids)
        plot_rd_gc(
            sub_win, sub_dp, rep_ids, genome_size, qc_dir,
            prefix, "window", "RD", ylim,
            gc_corr=gc_corr, gc_bin_median_std=gc_std,
            run_id=run_id, region_bed=region_bed, blacklist_bed=blacklist_bed,
            min_ylim=min_ylim,
        )
        return gc_corr, gc_std

    if dp_raw_ref is not None:
        _qc_plot_subset(dp_raw_ref, is_target_mask, "depth_raw_target")
        _qc_plot_subset(dp_raw_ref, ~is_target_mask, "depth_raw_antitarget")

    _qc_plot_subset(dp_corr_ref, is_target_mask, "depth_corrected_target")
    _qc_plot_subset(dp_corr_ref, ~is_target_mask, "depth_corrected_antitarget")

    pdf = PdfPages(stamp_path(os.path.join(qc_dir, "rd_correct.pdf"), run_id))
    if dp_raw_ref is not None:
        plot_gc_correction_pdf(gc_vals, dp_raw_ref, dp_corr_ref, rep_ids, pdf)
    else:
        plot_gc_correction_pdf(gc_vals, dp_corr_ref, dp_corr_ref, rep_ids, pdf)
    pdf.close()
else:
    logging.info("no GC data; skipping GC-related QC plots")
    ylim = max(np.nanquantile(dp_corr_ref, 0.99), 1.0) * 1.1
    min_ylim = min(np.nanquantile(dp_corr_ref, 0.01), 0.0) * 1.1
    plot_rd_gc(
        win_df, dp_corr_ref, rep_ids, genome_size, qc_dir,
        "depth_corrected", "window", "RD", ylim,
        run_id=run_id, region_bed=region_bed, blacklist_bed=blacklist_bed,
        min_ylim=min_ylim,
    )

del dp_raw_ref  # free memory before subdivision

win_df["_parent_chr"] = win_df["#CHR"].values
win_df["_parent_start"] = win_df["START"].values
win_df["_parent_end"] = win_df["END"].values

n_before_ext = len(win_df)
win_df = _extend_windows_to_midpoints(win_df, blacklist_bed, region_bed, genome_size)
total_gap_bp = int(
    (win_df["END"] - win_df["START"]).sum()
    - (win_df["_parent_end"] - win_df["_parent_start"]).sum()
)
logging.info(
    f"window extension: {n_before_ext} windows extended by {total_gap_bp:,} bp total"
)

if min_window_length > 0:
    sizes = (win_df["END"] - win_df["START"]).values
    n_pieces = np.where(
        sizes >= min_window_length, np.maximum(sizes // min_window_length, 1), 1
    ).astype(int)

    parent_idx = np.repeat(np.arange(len(win_df)), n_pieces)
    offsets = np.concatenate([np.arange(p) for p in n_pieces])
    orig_starts = win_df["START"].values[parent_idx]
    orig_ends = win_df["END"].values[parent_idx]

    new_starts = orig_starts + offsets * min_window_length
    is_last = np.zeros(len(parent_idx), dtype=bool)
    is_last[np.cumsum(n_pieces) - 1] = True
    new_ends = np.where(is_last, orig_ends, new_starts + min_window_length)

    win_df = win_df.iloc[parent_idx].reset_index(drop=True)
    win_df["START"] = new_starts
    win_df["END"] = new_ends
    logging.info(
        f"subdivision (min_window_length={min_window_length}): "
        f"{n_ref} → {len(win_df)} windows"
    )

if blacklist_bed is not None:
    n_before = len(win_df)
    pr_windows = pr.PyRanges(
        win_df.rename(columns={"#CHR": "Chromosome", "START": "Start", "END": "End"})
        .assign(_idx=np.arange(len(win_df)))
    )
    bl_idx = pr_windows.overlap(pr.read_bed(blacklist_bed)).df["_idx"].unique()
    keep_mask = ~np.isin(np.arange(len(win_df)), bl_idx)
    win_df = win_df.loc[keep_mask].reset_index(drop=True)
    logging.info(
        f"blacklist subtraction: {n_before} → {len(win_df)} windows "
        f"({n_before - len(win_df)} removed)"
    )

chrom_rank = {c: i for i, c in enumerate(CHROM_ORDER)}
win_df = win_df[win_df["#CHR"].isin(chrom_rank)].copy()
win_df["_chrom_rank"] = win_df["#CHR"].map(chrom_rank)
win_df = win_df.sort_values(["_chrom_rank", "START"]).reset_index(drop=True)
win_df = win_df.drop(columns=["_chrom_rank"])

n_windows = len(win_df)
logging.info(f"{n_windows} final windows after filtering and sorting")

parent_win = pd.DataFrame({
    "#CHR": win_df["_parent_chr"],
    "START": win_df["_parent_start"],
    "END": win_df["_parent_end"],
})
dp_corrected = _fill_depth_matrix(cnr_depth_maps, parent_win, fill=0.0)

win_df = win_df.drop(columns=["_parent_chr", "_parent_start", "_parent_end"])

n_zero = int((dp_corrected[:, 0] == 0.0).sum())
logging.info(
    f"depth filled: {n_windows - n_zero}/{n_windows} from .cnr, "
    f"{n_zero} filled with 0.0"
)
log_nan_summary("cnvkit corrected depth", dp_corrected, rep_ids, n_windows)

region_df = read_region_file(region_bed)
midpoints = ((win_df["START"] + win_df["END"]) // 2).values
chroms = win_df["#CHR"].values

r_chrs = region_df["#CHR"].values
r_starts = region_df["START"].values
r_ends = region_df["END"].values

region_ids = np.full(n_windows, None, dtype=object)
for k in range(len(region_df)):
    rid = f"{r_chrs[k]}:{r_starts[k]}-{r_ends[k]}"
    mask = (chroms == r_chrs[k]) & (midpoints >= r_starts[k]) & (midpoints < r_ends[k])
    region_ids[mask] = rid
win_df["region_id"] = region_ids
n_assigned = int((region_ids != None).sum())  # noqa: E711
logging.info(
    f"region_id assigned to {n_assigned}/{n_windows} "
    f"({n_assigned / max(n_windows, 1) * 100:.1f}%) windows"
)

np.savez_compressed(sm.output["dp_corrected"], mat=dp_corrected)

out_cols = ["#CHR", "START", "END", "region_id", "GC", "is_target"]
win_df[out_cols].to_csv(
    sm.output["window_df"],
    sep="\t",
    header=True,
    index=False,
    compression="gzip",
)

logging.info("finished cnvkit_to_window_dp.")
