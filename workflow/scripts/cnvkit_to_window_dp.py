"""Convert CNVkit .cnr files to window.dp.npz + window.tsv.gz.

Uses reference.cnn as the canonical window set (not .cnr, which drops windows
with bad coverage). Large windows are subdivided, blacklisted regions removed,
and target/antitarget annotated. Depth from .cnr is left-joined onto the
canonical windows, with unmatched windows filled with 0.0.

Raw coverage from .targetcoverage.cnn / .antitargetcoverage.cnn is loaded for
before/after comparison when available.
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
from io_utils import read_region_file

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages

setup_logging(sm.log[0])

# ---------------------------------------------------------------------------
# Natural chromosome sort order
# ---------------------------------------------------------------------------
CHROM_ORDER = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

# ---------------------------------------------------------------------------
# Inputs & params
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 1. Build canonical window set from reference.cnn
# ---------------------------------------------------------------------------
ref_cnn_file = sm.input["reference_cnn"]
ref_cnn = pd.read_table(ref_cnn_file, sep="\t")
ref_cnn = ref_cnn[ref_cnn["chromosome"].isin(target_chroms)].reset_index(drop=True)

# Assert no duplicate coordinates in reference.cnn
coords = ref_cnn[["chromosome", "start", "end"]].apply(tuple, axis=1)
assert coords.is_unique, "Duplicate coordinates in reference.cnn"

win_df = pd.DataFrame(
    {
        "#CHR": ref_cnn["chromosome"].values,
        "START": ref_cnn["start"].values,
        "END": ref_cnn["end"].values,
    }
)

# GC from reference.cnn
if "gc" in ref_cnn.columns:
    win_df["GC"] = ref_cnn["gc"].values
    logging.info(f"GC loaded from reference.cnn for {len(win_df)} windows")
else:
    win_df["GC"] = np.nan
    logging.warning(f"no gc column in reference.cnn ({ref_cnn_file})")

# is_target: gene not in ("Antitarget", "Background") → 1, else 0
if "gene" in ref_cnn.columns:
    win_df["is_target"] = (
        ~ref_cnn["gene"].isin(["Antitarget", "Background"])
    ).astype(np.int8)
else:
    win_df["is_target"] = np.int8(1)
    logging.warning("no gene column in reference.cnn; assuming all windows are target")

n_ref = len(win_df)
n_target = int(win_df["is_target"].sum())
logging.info(
    f"reference.cnn: {n_ref} windows ({n_target} target, {n_ref - n_target} antitarget)"
)

# ---------------------------------------------------------------------------
# 2. Load .cnr depth lookups (keyed by parent coordinates, before subdivision)
# ---------------------------------------------------------------------------
cnr_depth_maps = []
for rep_id in rep_ids:
    cnr_file = os.path.join(cnvkit_dir, f"{rep_id}.cnr")
    df = pd.read_table(cnr_file, sep="\t")
    df = df[df["chromosome"].isin(target_chroms)].reset_index(drop=True)
    depth_map = {}
    for _, row in df.iterrows():
        depth_map[(row["chromosome"], int(row["start"]), int(row["end"]))] = float(
            row["depth"]
        )
    cnr_depth_maps.append(depth_map)
    n_matched = sum(
        1
        for _, r in win_df.iterrows()
        if (r["#CHR"], int(r["START"]), int(r["END"])) in depth_map
    )
    logging.info(
        f"  {rep_id}: {len(df)} .cnr bins, {n_matched}/{n_ref} matched to reference.cnn"
    )

# Store parent coordinates for depth lookup after subdivision
win_df["_parent_chr"] = win_df["#CHR"].values
win_df["_parent_start"] = win_df["START"].values
win_df["_parent_end"] = win_df["END"].values

# ---------------------------------------------------------------------------
# 3. Subdivide large windows
# ---------------------------------------------------------------------------
if min_window_length > 0:
    rows = []
    for _, row in win_df.iterrows():
        size = int(row["END"] - row["START"])
        if size >= min_window_length:
            n_pieces = size // min_window_length
            if n_pieces < 1:
                n_pieces = 1
            start = int(row["START"])
            for j in range(n_pieces):
                sub_start = start + j * min_window_length
                if j == n_pieces - 1:
                    sub_end = int(row["END"])
                else:
                    sub_end = sub_start + min_window_length
                new_row = row.copy()
                new_row["START"] = sub_start
                new_row["END"] = sub_end
                rows.append(new_row)
        else:
            rows.append(row)

    win_df_sub = pd.DataFrame(rows).reset_index(drop=True)
    logging.info(
        f"subdivision (min_window_length={min_window_length}): "
        f"{n_ref} → {len(win_df_sub)} windows"
    )
    win_df = win_df_sub

# ---------------------------------------------------------------------------
# 4. Subtract blacklisted regions
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 5. Sort by natural chromosome order, then START
# ---------------------------------------------------------------------------
chrom_rank = {c: i for i, c in enumerate(CHROM_ORDER)}
win_df = win_df[win_df["#CHR"].isin(chrom_rank)].copy()
win_df["_chrom_rank"] = win_df["#CHR"].map(chrom_rank)
win_df = win_df.sort_values(["_chrom_rank", "START"]).reset_index(drop=True)
win_df = win_df.drop(columns=["_chrom_rank"])

n_windows = len(win_df)
logging.info(f"{n_windows} final windows after filtering and sorting")

# ---------------------------------------------------------------------------
# 6. Build dp_corrected matrix (depth from parent .cnr, or 0.0 if missing)
# ---------------------------------------------------------------------------
dp_corrected = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, depth_map in enumerate(cnr_depth_maps):
    for j in range(n_windows):
        key = (
            win_df.iloc[j]["_parent_chr"],
            int(win_df.iloc[j]["_parent_start"]),
            int(win_df.iloc[j]["_parent_end"]),
        )
        dp_corrected[j, i] = depth_map.get(key, 0.0)

n_zero = int((dp_corrected[:, 0] == 0.0).sum())
logging.info(
    f"depth filled: {n_windows - n_zero}/{n_windows} from .cnr, "
    f"{n_zero} filled with 0.0"
)

# Drop parent coordinate columns
win_df = win_df.drop(columns=["_parent_chr", "_parent_start", "_parent_end"])

log_nan_summary("cnvkit corrected depth", dp_corrected, rep_ids, n_windows)

# ---------------------------------------------------------------------------
# 7. Region ID assignment
# ---------------------------------------------------------------------------
if region_bed is not None:
    region_df = read_region_file(region_bed)
    midpoints = ((win_df["START"] + win_df["END"]) // 2).values
    chroms = win_df["#CHR"].values

    region_ids = np.full(n_windows, None, dtype=object)
    for _, row in region_df.iterrows():
        r_chr = row["#CHR"] if "#CHR" in region_df.columns else row["Chromosome"]
        r_start = row["START"] if "START" in region_df.columns else row["Start"]
        r_end = row["END"] if "END" in region_df.columns else row["End"]
        rid = f"{r_chr}:{r_start}-{r_end}"
        mask = (chroms == r_chr) & (midpoints >= r_start) & (midpoints < r_end)
        region_ids[mask] = rid
    win_df["region_id"] = region_ids
    n_assigned = int((region_ids != None).sum())  # noqa: E711
    logging.info(
        f"region_id assigned to {n_assigned}/{n_windows} "
        f"({n_assigned / max(n_windows, 1) * 100:.1f}%) windows"
    )
else:
    win_df["region_id"] = 0
    logging.info("no region_bed; setting region_id=0 for all windows")

# ---------------------------------------------------------------------------
# 8. Load raw coverage for QC plots (optional)
# ---------------------------------------------------------------------------
dp_raw = None
cnr_coords = win_df[["#CHR", "START", "END"]]
try:
    dp_raw = np.full((n_windows, nsamples), np.nan, dtype=np.float32)
    for i, rep_id in enumerate(rep_ids):
        tgt_file = os.path.join(cnvkit_dir, f"{rep_id}.targetcoverage.cnn")
        anti_file = os.path.join(cnvkit_dir, f"{rep_id}.antitargetcoverage.cnn")
        tgt_df = pd.read_table(tgt_file, sep="\t")
        anti_df = pd.read_table(anti_file, sep="\t")
        raw_df = pd.concat([tgt_df, anti_df], ignore_index=True)
        raw_df = raw_df[raw_df["chromosome"].isin(target_chroms)].reset_index(drop=True)
        raw_df = raw_df.rename(columns={"chromosome": "#CHR", "start": "START", "end": "END"})
        merged = cnr_coords.merge(
            raw_df[["#CHR", "START", "END", "depth"]],
            on=["#CHR", "START", "END"],
            how="left",
        )
        dp_raw[:, i] = merged["depth"].to_numpy(dtype=np.float32)
    n_matched = int(np.isfinite(dp_raw[:, 0]).sum())
    logging.info(
        f"loaded raw coverage from .cnn files: {n_matched}/{n_windows} "
        f"({n_matched / max(n_windows, 1) * 100:.1f}%) bins matched"
    )
    if n_matched == 0:
        dp_raw = None
except Exception as e:
    logging.warning(
        f"could not load raw .cnn files: {e}; skipping before/after comparison"
    )

# ---------------------------------------------------------------------------
# 9. QC plots
# ---------------------------------------------------------------------------
gc_vals = win_df["GC"].to_numpy()
has_gc = np.isfinite(gc_vals).any()

if has_gc:
    if dp_raw is not None:
        rd_raw_ylim = max(np.nanquantile(dp_raw, 0.99), 1.0) * 1.1
        gc_corr_before, gc_std_before = compute_gc_rd_stats(dp_raw, gc_vals, rep_ids)
        plot_rd_gc(
            win_df,
            dp_raw,
            rep_ids,
            genome_size,
            qc_dir,
            "depth_before_correction",
            "window",
            "RD",
            rd_raw_ylim,
            gc_corr=gc_corr_before,
            gc_bin_median_std=gc_std_before,
            run_id=run_id,
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
        )

    rd_ylim = max(np.nanquantile(dp_corrected, 0.99), 1.0) * 1.1
    gc_corr_after, gc_std_after = compute_gc_rd_stats(dp_corrected, gc_vals, rep_ids)
    plot_rd_gc(
        win_df,
        dp_corrected,
        rep_ids,
        genome_size,
        qc_dir,
        "depth_cnvkit_corrected",
        "window",
        "RD",
        rd_ylim,
        gc_corr=gc_corr_after,
        gc_bin_median_std=gc_std_after,
        run_id=run_id,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
    )

    pdf = PdfPages(stamp_path(os.path.join(qc_dir, "rd_correct.pdf"), run_id))
    if dp_raw is not None:
        plot_gc_correction_pdf(gc_vals, dp_raw, dp_corrected, rep_ids, pdf)
    else:
        plot_gc_correction_pdf(gc_vals, dp_corrected, dp_corrected, rep_ids, pdf)
    pdf.close()
else:
    logging.info("no GC data available; skipping GC-related QC plots")
    rd_ylim = max(np.nanquantile(dp_corrected, 0.99), 1.0) * 1.1
    plot_rd_gc(
        win_df,
        dp_corrected,
        rep_ids,
        genome_size,
        qc_dir,
        "depth_cnvkit_corrected",
        "window",
        "RD",
        rd_ylim,
        run_id=run_id,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
    )

# ---------------------------------------------------------------------------
# 10. Write outputs (no NaN filtering)
# ---------------------------------------------------------------------------
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
