"""Per-window LOWESS bias correction for bulk samples.

Loads mosdepth fixed-window depth, joins with pre-filtered window BED
(GC/MAP/REPLI/region_id/is_target), applies correct_readcount() per
sample, and saves corrected depth matrix + filtered window dataframe.

The window BED is expected to be pre-filtered by region and blacklist
(produced by build_window_bed.py), with region_id and optional is_target
columns already present.
"""

import os, logging
from snakemake.script import snakemake as sm

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from utils import setup_logging
from count_reads_utils import (
    log_nan_summary,
    log_mad_and_plot,
)
from rd_correct_utils import (
    correct_readcount,
    correct_readcount_wes,
    plot_gc_correction_pdf,
)

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages

setup_logging(sm.log[0])

# ---------------------------------------------------------------------------
# Inputs / params
# ---------------------------------------------------------------------------
sample_file = sm.input["sample_file"]
window_bed_file = sm.input["window_bed"]
genome_size = sm.input["genome_size"]
assay_type = str(sm.params["assay_type"])

mosdepth_dir = sm.params["mosdepth_dir"]
chromosomes = sm.params["chromosomes"]
samplesize = int(sm.params["samplesize"])
routlier = float(sm.params["routlier"])
doutlier = float(sm.params["doutlier"])
min_mappability = float(sm.params["min_mappability"])
gc_correct = bool(sm.params["gc_correct"])
rt_correct = bool(sm.params["rt_correct"])

qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Load sample info
# ---------------------------------------------------------------------------
sample_df = pd.read_table(sample_file, sep="\t")
rep_ids = sample_df["REP_ID"].astype(str).tolist()
nsamples = len(sample_df)
target_chroms = {f"chr{c}" for c in chromosomes}
join_keys = ["#CHR", "START", "END"]

logging.info(f"rd_correct: {nsamples} samples, {len(target_chroms)} chroms")

# ---------------------------------------------------------------------------
# Load window BED (pre-filtered by region and blacklist)
# ---------------------------------------------------------------------------
logging.info("load window BED and mosdepth depth")
win_df = pd.read_table(window_bed_file, sep="\t")
assert "#CHR" in win_df.columns and "GC" in win_df.columns, (
    f"window_bed must have #CHR, START, END, GC columns; got {win_df.columns.tolist()}"
)

# Filter to target chromosomes
win_df = win_df[win_df["#CHR"].isin(target_chroms)].reset_index(drop=True)

# ---------------------------------------------------------------------------
# Load mosdepth depth and join with window BED
# ---------------------------------------------------------------------------
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

# Merge mosdepth coords with window BED to get covariates
win_df = pd.merge(left=coords, right=win_df, on=join_keys, how="left", sort=False)
_gc_matched = int(win_df["GC"].notna().sum())
logging.info(
    f"GC BED matched {_gc_matched}/{n_windows} ({_gc_matched / max(n_windows, 1) * 100:.1f}%)"
)

dp_raw = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, mos_df in enumerate(mos_dfs):
    dp_raw[:, i] = mos_df["DEPTH"].to_numpy(dtype=np.float32)

rd_raw_ylim = max(np.nanquantile(dp_raw, 0.99), 1.0)
log_mad_and_plot(
    win_df,
    dp_raw,
    rep_ids,
    genome_size,
    qc_dir,
    "depth_before_correction",
    "window",
    "RD",
    rd_raw_ylim,
)

logging.info(f"{n_windows} windows for bias correction")

# ---------------------------------------------------------------------------
# Use is_target from window BED (WES) or default to all-target (WGS)
# ---------------------------------------------------------------------------
if "is_target" in win_df.columns:
    is_target = win_df["is_target"].astype(bool).to_numpy()
    n_tgt = int(is_target.sum())
    n_anti = n_windows - n_tgt
    logging.info(
        f"WES window labels from window BED: {n_tgt} target, {n_anti} antitarget"
    )
else:
    is_target = np.ones(n_windows, dtype=bool)
    win_df["is_target"] = True

# ---------------------------------------------------------------------------
# Extract bias covariates
# ---------------------------------------------------------------------------
gc_vals = win_df["GC"].to_numpy()
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

# ---------------------------------------------------------------------------
# Bias correction (separate for target / antitarget when WES)
# ---------------------------------------------------------------------------
if gc_correct:
    dp_corrected = np.zeros_like(dp_raw, dtype=np.float32)

    if assay_type == "bulkWES":
        logging.info("applying correct_readcount_wes per sample")
        for i, rep_id in enumerate(rep_ids):
            logging.info(f"  correcting {rep_id}")
            dp_corrected[:, i] = correct_readcount_wes(
                dp_raw[:, i],
                gc_vals,
                is_target,
            )
    else:
        logging.info("applying correct_readcount per sample")
        for i, rep_id in enumerate(rep_ids):
            logging.info(f"correcting {rep_id}")
            dp_corrected[:, i] = correct_readcount(
                dp_raw[:, i],
                gc_vals,
                mappability=map_vals,
                repliseq=repli_vals,
                samplesize=samplesize,
                routlier=routlier,
                doutlier=doutlier,
                min_mappability=min_mappability,
            )
else:
    logging.info("gc_correct=False; skipping bias correction")
    dp_corrected = dp_raw.copy()

log_nan_summary("corrected depth", dp_corrected, rep_ids, n_windows)

pdf = PdfPages(os.path.join(qc_dir, "depth_correction.pdf"))
plot_gc_correction_pdf(gc_vals, dp_raw, dp_corrected, rep_ids, pdf)
pdf.close()

rd_ylim = max(np.nanquantile(dp_corrected, 0.99), 1.0)
log_mad_and_plot(
    win_df,
    dp_corrected,
    rep_ids,
    genome_size,
    qc_dir,
    "depth_after_correction",
    "window",
    "RD",
    rd_ylim,
    smooth=True,
)

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
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
