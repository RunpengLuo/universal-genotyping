"""Per-window LOWESS bias correction for bulk samples.

Loads mosdepth fixed-window depth, joins with pre-filtered window BED
(GC/MAP/REPLI/region_id), applies correct_readcount_lowess() per
sample, and saves corrected depth matrix + filtered window dataframe.

The window BED is expected to be pre-filtered by region and blacklist
(produced by build_window_bed.py), with region_id column already present.
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
from io_utils import compute_depth_statistics
from count_reads_utils import compute_gc_rd_stats
from rd_correct_utils import (
    correct_readcount_lowess,
    correct_readcount_quadreg,
)
from plot_utils import plot_rd_gc, plot_gc_correction_pdf

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages

setup_logging(sm.log[0])

sample_file = sm.input["sample_file"]
window_bed_file = sm.input["window_bed"]
genome_size = sm.input["genome_size"]
region_bed = sm.input["region_bed"] or None
blacklist_bed = maybe_path(sm.input.get("blacklist_bed", None))
assay_type = str(sm.params["assay_type"])

mosdepth_dir = sm.params["mosdepth_dir"]
chromosomes = sm.params["chromosomes"]
samplesize = int(sm.params["samplesize"])
routlier = float(sm.params["routlier"])
doutlier = float(sm.params["doutlier"])
min_mappability = float(sm.params["min_mappability"])
gc_correct = bool(sm.params["gc_correct"])
gc_correct_method = str(sm.params.get("gc_correct_method", "median"))
rt_correct = bool(sm.params["rt_correct"])
qc_dir = sm.params["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
run_id = getattr(sm.params, "run_id", "")

sample_df = pd.read_table(sample_file, sep="\t")
rep_ids = sample_df["REP_ID"].astype(str).tolist()
nsamples = len(sample_df)
target_chroms = {f"chr{c}" for c in chromosomes}
join_keys = ["#CHR", "START", "END"]

logging.info(f"rd_correct: {nsamples} samples, {len(target_chroms)} chroms")

logging.info("load window BED and mosdepth depth")
win_df = pd.read_table(window_bed_file, sep="\t")
assert "#CHR" in win_df.columns and "GC" in win_df.columns, (
    f"window_bed must have #CHR, START, END, GC columns; got {win_df.columns.tolist()}"
)

win_df = win_df[win_df["#CHR"].isin(target_chroms)].reset_index(drop=True)

mos_dfs = []
for rep_id in rep_ids:
    mos_file = os.path.join(mosdepth_dir, f"{rep_id}.regions.bed.gz")
    mos_df = pd.read_table(
        mos_file, sep="\t", header=None, names=["#CHR", "START", "END", "DEPTH"]
    )
    mos_df = mos_df[mos_df["#CHR"].isin(target_chroms)].reset_index(drop=True)
    mos_dfs.append(mos_df)

coords = mos_dfs[0][join_keys].copy()
n_windows = len(coords)
logging.info(f"{n_windows} windows across {len(target_chroms)} chromosomes")

win_df = pd.merge(left=coords, right=win_df, on=join_keys, how="left", sort=False)
_gc_matched = int(win_df["GC"].notna().sum())
logging.info(
    f"GC BED matched {_gc_matched}/{n_windows} ({_gc_matched / max(n_windows, 1) * 100:.1f}%)"
)

dp_raw = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, mos_df in enumerate(mos_dfs):
    dp_raw[:, i] = mos_df["DEPTH"].to_numpy(dtype=np.float32)

sample_ids = sample_df["SAMPLE"].astype(str).tolist()
depth_stats = compute_depth_statistics(dp_raw, win_df, sample_ids)
depth_stats.to_csv(sm.output["depth_stats"], sep="\t", index=False)
logging.info(f"wrote depth statistics to {sm.output['depth_stats']}")
for _, row in depth_stats[depth_stats["#CHR"] == "TOTAL"].iterrows():
    logging.info(
        f"  {row['SAMPLE']}: mean={row['mean_depth']:.2f}, median={row['median_depth']:.2f}"
    )

snps = pd.read_table(sm.input["snp_info"], sep="\t")
tot_mtx_snp = np.load(sm.input["tot_mtx_snp"])["mat"].astype(np.int32)
a_mtx_snp = np.load(sm.input["a_mtx_snp"])["mat"].astype(np.int32)
b_mtx_snp = np.load(sm.input["b_mtx_snp"])["mat"].astype(np.int32)
logging.info(f"loaded {len(snps)} SNPs, allele matrices shape={tot_mtx_snp.shape}")

gc_vals = win_df["GC"].to_numpy()

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

logging.info(f"{n_windows} windows for bias correction")

map_vals = win_df["MAP"].to_numpy() if gc_correct and "MAP" in win_df.columns else None
repli_vals = (
    win_df["REPLI"].to_numpy(dtype=np.float64)
    if rt_correct and "REPLI" in win_df.columns
    else None
)
if repli_vals is not None:
    _n_finite = int(np.isfinite(repli_vals).sum())
    logging.info(
        f"REPLI column: {_n_finite}/{n_windows} ({_n_finite / max(n_windows, 1) * 100:.1f}%) finite"
    )
else:
    logging.info("no REPLI column; skipping replication timing correction")

gc_rmse_list = None
if gc_correct:
    dp_corrected = np.zeros_like(dp_raw, dtype=np.float32)
    gc_rmse_list = []

    if gc_correct_method == "median":
        logging.info("applying correct_readcount_quadreg per sample")
        for i, rep_id in enumerate(rep_ids):
            logging.info(f"correcting {rep_id}")
            dp_corrected[:, i], gc_rmse = correct_readcount_quadreg(
                dp_raw[:, i],
                gc_vals,
                mappability=map_vals,
                repliseq=repli_vals,
                doutlier=doutlier,
                min_mappability=min_mappability,
            )
            gc_rmse_list.append(gc_rmse)
    else:
        logging.info("applying correct_readcount_lowess per sample")
        for i, rep_id in enumerate(rep_ids):
            logging.info(f"correcting {rep_id}")
            dp_corrected[:, i], gc_rmse = correct_readcount_lowess(
                dp_raw[:, i],
                gc_vals,
                mappability=map_vals,
                repliseq=repli_vals,
                samplesize=samplesize,
                routlier=routlier,
                doutlier=doutlier,
                min_mappability=min_mappability,
            )
            gc_rmse_list.append(gc_rmse)
else:
    logging.info("gc_correct=False; skipping bias correction")
    dp_corrected = dp_raw.copy()

pdf = PdfPages(stamp_path(os.path.join(qc_dir, "rd_correct.pdf"), run_id))
plot_gc_correction_pdf(
    gc_vals,
    dp_raw,
    dp_corrected,
    rep_ids,
    pdf,
    gc_rmse=gc_rmse_list,
    mappability=map_vals,
    repliseq=repli_vals,
)
pdf.close()

rd_ylim = max(np.nanquantile(dp_corrected, 0.99), 1.0) * 1.1
gc_corr_after, gc_std_after = compute_gc_rd_stats(dp_corrected, gc_vals, rep_ids)
plot_rd_gc(
    win_df,
    dp_corrected,
    rep_ids,
    genome_size,
    qc_dir,
    "depth_after_correction",
    "window",
    "RD",
    rd_ylim,
    gc_corr=gc_corr_after,
    gc_bin_median_std=gc_std_after,
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
if "MAP" in win_df.columns:
    out_cols.append("MAP")
if "REPLI" in win_df.columns:
    out_cols.append("REPLI")
win_df[out_cols].to_csv(
    sm.output["window_df"],
    sep="\t",
    header=True,
    index=False,
    compression="gzip",
)

logging.info("finished rd_correct.")
