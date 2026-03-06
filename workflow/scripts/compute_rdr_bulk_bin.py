import os, sys, gzip, logging
from snakemake.script import snakemake as sm

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from utils import *
from bias_correction import (
    gc_correct_depth_loess,
    plot_depth_gc_correction,
    bias_correction_rdr_quantreg,
    bias_correction_rdr_spline,
)
from postprocess_utils import (
    plot_1d_sample,
    log_nan_summary,
    log_mad_and_plot,
    compute_overlap_weights,
    weighted_bincount_mean,
)

##################################################
logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# input
sample_file = sm.input["sample_file"]
bb_file = sm.input["bb_file"]
bias_bed_file = sm.input["bias_bed"]
genome_size = sm.input["genome_size"]
region_bed_file = maybe_path(sm.input["region_bed"])

baf_mtx_bb = np.load(sm.input["baf_mtx_bb"])["mat"]

sample_name = sm.params["sample_name"]
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
mosdepth_dir = sm.params["mosdepth_dir"]
correction_method = str(sm.params["correction_method"])
gc_correct = bool(sm.params["gc_correct"])
rt_correct = bool(sm.params["rt_correct"])
chromosomes = sm.params["chromosomes"]

logging.info(f"run compute_rdr_bulk (gc_correct={gc_correct}, rt_correct={rt_correct})")

##################################################
logging.info("load bias BED")
bias_bed_full = pd.read_table(bias_bed_file, sep="\t")
assert "#CHR" in bias_bed_full.columns and "GC" in bias_bed_full.columns, (
    f"bias_bed must have #CHR, START, END, GC columns; got {bias_bed_full.columns.tolist()}"
)

##################################################
logging.info("concat mosdepth per-sample depth data")
sample_df = pd.read_table(sample_file, sep="\t")
rep_ids = sample_df["REP_ID"].astype(str).tolist()
sample_types = sample_df["sample_type"].tolist()
bbs = pd.read_table(bb_file, sep="\t")
nsamples = len(sample_df)
ntumor_samples = (sample_df["sample_type"] == "tumor").sum()
nbb = len(bbs)
for rep_id in rep_ids:
    mos_file = os.path.join(mosdepth_dir, f"{rep_id}.regions.bed.gz")
    mos_df = pd.read_table(
        mos_file, sep="\t", header=None, names=["#CHR", "START", "END", rep_id]
    )
    dup = mos_df.duplicated(subset=["#CHR", "START", "END"], keep=False)
    assert not dup.any(), f"{mos_file} has duplicate intervals: {int(dup.sum())}"
    assert len(mos_df) == len(bbs), (
        f"{mos_file} rowcount={len(mos_df)} != bbs rowcount={len(bbs)}"
    )
    bbs = bbs.merge(mos_df, how="left", on=["#CHR", "START", "END"], sort=False)
dp_mtx_bb = bbs[rep_ids].to_numpy(dtype=np.float32)
np.savez_compressed(sm.output["dp_mtx_bb"], mat=dp_mtx_bb)

##################################################
# Aggregate bias_bed windows into adaptive bins via overlap weighting
logging.info("aggregate bias_bed windows into adaptive bins")
target_chroms = {f"chr{c}" for c in chromosomes}
bias_bed_full = bias_bed_full[bias_bed_full["#CHR"].isin(target_chroms)].reset_index(drop=True)

bias_cols = ["GC"]
if "MAP" in bias_bed_full.columns:
    bias_cols.append("MAP")
if "REPLI" in bias_bed_full.columns:
    bias_cols.append("REPLI")

bias_arrays = {col: bias_bed_full[col].to_numpy(dtype=np.float64) for col in bias_cols}
bb_bias = {col: np.full(nbb, np.nan, dtype=np.float64) for col in bias_cols}

win_chr = bias_bed_full["#CHR"].to_numpy()
win_start = bias_bed_full["START"].to_numpy(dtype=np.int64)
win_end = bias_bed_full["END"].to_numpy(dtype=np.int64)

for chrom, bb_grp in bbs.groupby("#CHR", sort=False):
    win_mask = win_chr == chrom
    win_idx_chr = np.where(win_mask)[0]
    if len(win_idx_chr) == 0:
        continue

    bb_starts = bb_grp["START"].to_numpy(dtype=np.int64)
    bb_ends = bb_grp["END"].to_numpy(dtype=np.int64)
    bb_global_idx = bb_grp.index.to_numpy()
    n_bb_chr = len(bb_grp)

    ov_win, ov_bin, ov_wt = compute_overlap_weights(
        win_start[win_idx_chr],
        win_end[win_idx_chr],
        bb_starts,
        bb_ends,
    )
    if len(ov_win) == 0:
        logging.info(f"  {chrom}: {n_bb_chr} bins, 0 windows assigned")
        continue

    for col in bias_cols:
        vals = bias_arrays[col][win_idx_chr][ov_win]
        bb_bias[col][bb_global_idx] = weighted_bincount_mean(
            vals, ov_bin, ov_wt, n_bb_chr
        )

    logging.info(f"  {chrom}: {n_bb_chr} bins, {len(np.unique(ov_win))} windows assigned")

corr_factors = bbs[["#CHR", "START", "END"]].copy()
for col in bias_cols:
    corr_factors[col] = bb_bias[col]
corr_factors.to_csv(sm.output["corr_factors"], sep="\t", header=True, index=False)
logging.info(f"wrote corr_factors: {bias_cols}, {nbb} rows")

##################################################
has_normal = "normal" in sample_types
logging.info(f"compute RDRs, has_normal={has_normal}")
assert has_normal, "no normal sample, TODO"
tumor_sidx = {False: 0, True: 1}[has_normal]

gc_vals = corr_factors["GC"].to_numpy()
mapp_vals = corr_factors["MAP"].to_numpy() if gc_correct and "MAP" in corr_factors.columns else None
repli_vals = (
    corr_factors["REPLI"].to_numpy(dtype=np.float64)
    if rt_correct and "REPLI" in corr_factors.columns
    else None
)

rd_raw_ylim = max(np.nanquantile(dp_mtx_bb, 0.99), 1.0)
log_nan_summary("raw depth", dp_mtx_bb, rep_ids, nbb)
log_mad_and_plot(
    bbs, dp_mtx_bb, rep_ids, genome_size, qc_dir,
    "depth_before_correction", "bb", "RD", rd_raw_ylim,
    smooth=True, region_bed=region_bed_file,
)

##################################################
# Loess: correct raw depth for GC bias BEFORE computing RDR
if gc_correct and correction_method == "loess":
    logging.info("loess GC correction on raw depth (pre-RDR)")

    dp_mtx_bb_raw = dp_mtx_bb.copy()  # keep raw for diagnostics
    dp_mtx_bb = gc_correct_depth_loess(
        dp_mtx_bb,
        gc_vals,
        mappability=mapp_vals,
    )

    # diagnostic plots: depth vs GC before/after for each sample
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages(os.path.join(qc_dir, "loess_depth_correction.pdf"))
    for i, rep_id in enumerate(rep_ids):
        plot_depth_gc_correction(
            gc_vals,
            dp_mtx_bb_raw[:, i],
            dp_mtx_bb[:, i],
            rep_id,
            pdf,
        )
    pdf.close()

    log_nan_summary("corrected depth", dp_mtx_bb, rep_ids, nbb)
    rd_corr_ylim = max(np.nanquantile(dp_mtx_bb, 0.99), 1.0)
    log_mad_and_plot(
        bbs, dp_mtx_bb, rep_ids, genome_size, qc_dir,
        "depth_after_correction", "bb", "RD", rd_corr_ylim,
        smooth=True, region_bed=region_bed_file,
    )

##################################################
if has_normal:
    bases_mtx_bb = dp_mtx_bb * bbs["BLOCKSIZE"].to_numpy()[:, None]
    total_bases = np.sum(bases_mtx_bb, axis=0)
    library_correction = total_bases[0] / total_bases[1:]
    logging.info(f"RDR library normalization factor: {library_correction}")

    # view
    normal_dp = dp_mtx_bb[:, 0]

    mask = np.isfinite(normal_dp) & (normal_dp > 0)
    num_valid = int(mask.sum())

    if num_valid < nbb:
        logging.warning(
            f"#invalid normal sample bins due to infinite/zero depth={nbb - num_valid}/{nbb}"
        )
        if num_valid == 0:
            raise ValueError("All normal bins invalid; cannot impute.")
        fill = float(np.nanmedian(normal_dp[mask]))
        logging.warning(
            f"impute normal depth via global nanmedian over valid normal bins: {fill}"
        )
        normal_dp[~mask] = fill

    rdr_mtx_bb = dp_mtx_bb[:, 1:] / dp_mtx_bb[:, 0][:, None]
    rdr_mtx_bb *= library_correction[None, :]
else:
    # TODO panel of normal PON, metin longread branch.
    raise ValueError("no normal sample, TODO")

tumor_rep_ids = rep_ids[tumor_sidx:]
rdr_ylim = np.round(np.nanquantile(rdr_mtx_bb, 0.99)).astype(int) + 1
log_nan_summary("RDR (pre-corr)", rdr_mtx_bb, tumor_rep_ids, nbb)
log_mad_and_plot(
    bbs, rdr_mtx_bb, tumor_rep_ids, genome_size, qc_dir,
    "rdr_before_correction", "bb", "RDR", rdr_ylim,
    smooth=True, region_bed=region_bed_file,
)

##################################################
# Post-RDR GC correction (quantile / spline methods — skip if loess already applied)
if gc_correct and correction_method != "loess":
    has_mapp = mapp_vals is not None

    if correction_method == "spline":
        rdr_mtx_bb = bias_correction_rdr_spline(
            rdr_mtx_bb,
            corr_factors,
            rep_ids[tumor_sidx:],
            has_mapp=has_mapp,
            rt_vals=repli_vals,
            out_dir=qc_dir,
        )
    else:  # "quantile" — current default
        rdr_mtx_bb = bias_correction_rdr_quantreg(
            rdr_mtx_bb,
            corr_factors,
            rep_ids[tumor_sidx:],
            has_mapp=has_mapp,
            rt_vals=repli_vals,
            out_dir=qc_dir,
        )

# bias covariate plots
for col in bias_cols:
    plot_1d_sample(
        corr_factors,
        corr_factors[col].to_numpy(),
        genome_size,
        os.path.join(qc_dir, f"{col}.pdf"),
        unit="bb",
        val_type=col,
        smooth=True,
        region_bed=region_bed_file,
    )

# plot per-sample RDRs (after corrections)
log_nan_summary("RDR (post-corr)", rdr_mtx_bb, tumor_rep_ids, nbb)
log_mad_and_plot(
    bbs, rdr_mtx_bb, tumor_rep_ids, genome_size, qc_dir,
    "rdr", "bb", "RDR", rdr_ylim,
    smooth=True, region_bed=region_bed_file,
)

np.savez_compressed(sm.output["rdr_mtx_bb"], mat=rdr_mtx_bb)
logging.info(f"finished.")
