"""Convert CNVkit .cnr files to window.dp.npz + window.tsv.gz.

Reads per-replicate .cnr files (from cnvkit fix), builds a unified window
DataFrame with region_id assignments, and generates QC plots. Raw coverage
from .targetcoverage.cnn / .antitargetcoverage.cnn is loaded for before/after
comparison when available.
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

sample_file = sm.input["sample_file"]
region_bed = sm.input["region_bed"] or None
blacklist_bed = maybe_path(sm.input.get("blacklist_bed", None))
genome_size = sm.input["genome_size"]

cnvkit_dir = sm.params["cnvkit_dir"]
chromosomes = sm.params["chromosomes"]

qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
run_id = getattr(sm.params, "run_id", "")

sample_df = pd.read_table(sample_file, sep="\t")
rep_ids = sample_df["REP_ID"].astype(str).tolist()
nsamples = len(sample_df)
target_chroms = {f"chr{c}" for c in chromosomes}

logging.info(f"cnvkit_to_window_dp: {nsamples} samples, {len(target_chroms)} chroms")

cnr_dfs = []
for rep_id in rep_ids:
    cnr_file = os.path.join(cnvkit_dir, f"{rep_id}.cnr")
    df = pd.read_table(cnr_file, sep="\t")
    df = df[df["chromosome"].isin(target_chroms)].reset_index(drop=True)
    cnr_dfs.append(df)
    logging.info(f"  loaded {len(df)} bins from {cnr_file}")

ref_cnr = cnr_dfs[0]
n_windows = len(ref_cnr)
logging.info(f"{n_windows} windows across {len(target_chroms)} chromosomes")

for i, df in enumerate(cnr_dfs[1:], 1):
    assert len(df) == n_windows, (
        f"sample {rep_ids[i]} has {len(df)} bins vs {n_windows} in {rep_ids[0]}"
    )
    assert (df["chromosome"].values == ref_cnr["chromosome"].values).all(), (
        f"chromosome mismatch between {rep_ids[i]} and {rep_ids[0]}"
    )
    assert (df["start"].values == ref_cnr["start"].values).all(), (
        f"start coordinate mismatch between {rep_ids[i]} and {rep_ids[0]}"
    )

win_df = pd.DataFrame(
    {
        "#CHR": ref_cnr["chromosome"].values,
        "START": ref_cnr["start"].values,
        "END": ref_cnr["end"].values,
    }
)

if "gc" in ref_cnr.columns:
    win_df["GC"] = ref_cnr["gc"].values
    logging.info("GC column found in .cnr")
else:
    # Load GC from reference.cnn (cnvkit fix drops gc/rmask columns from .cnr)
    ref_cnn_file = sm.input["reference_cnn"]
    ref_cnn = pd.read_table(ref_cnn_file, sep="\t")
    ref_cnn = ref_cnn[ref_cnn["chromosome"].isin(target_chroms)].reset_index(drop=True)
    if "gc" in ref_cnn.columns:
        merged = win_df[["#CHR", "START", "END"]].merge(
            ref_cnn.rename(columns={"chromosome": "#CHR", "start": "START", "end": "END"})[
                ["#CHR", "START", "END", "gc"]
            ],
            on=["#CHR", "START", "END"],
            how="left",
        )
        win_df["GC"] = merged["gc"].values
        n_matched = int(merged["gc"].notna().sum())
        logging.info(
            f"GC loaded from reference.cnn: {n_matched}/{n_windows} "
            f"({n_matched / max(n_windows, 1) * 100:.1f}%) bins matched"
        )
    else:
        win_df["GC"] = np.nan
        logging.warning(f"no gc column in reference.cnn ({ref_cnn_file})")

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

dp_corrected = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, df in enumerate(cnr_dfs):
    dp_corrected[:, i] = df["depth"].to_numpy(dtype=np.float32)

log_nan_summary("cnvkit corrected depth", dp_corrected, rep_ids, n_windows)

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

nan_mask = np.isnan(dp_corrected).any(axis=1)
n_nan_rows = int(nan_mask.sum())
n_valid = n_windows - n_nan_rows
logging.info(
    f"NaN row filter: {n_nan_rows}/{n_windows} windows have NaN, "
    f"keeping {n_valid} ({n_valid / max(n_windows, 1) * 100:.1f}%)"
)

if n_nan_rows > 0:
    valid = ~nan_mask
    dp_corrected = dp_corrected[valid]
    win_df = win_df.loc[valid].reset_index(drop=True)

np.savez_compressed(sm.output["dp_corrected"], mat=dp_corrected)

out_cols = ["#CHR", "START", "END", "region_id", "GC"]
win_df[out_cols].to_csv(
    sm.output["window_df"],
    sep="\t",
    header=True,
    index=False,
    compression="gzip",
)

logging.info("finished cnvkit_to_window_dp.")
