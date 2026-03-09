"""SNP-informed adaptive binning + window→bin depth aggregation + RDR.

Combines:
1. Adaptive binning from SNP data (switch probs, MSR/MSPB)
2. Aggregation of corrected window depth into adaptive bins
3. RDR computation (tumor/normal ratio or median-centering)
"""

import os
import shutil
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
from aggregation_utils import (
    adaptive_binning_windows,
    assign_pos_to_range,
    matrix_segmentation,
)
from combine_counts_utils import plot_allele_freqs
from count_reads_utils import log_nan_summary, plot_1d_multi_sample
from switchprobs import (
    interp_cM_blocks,
    estimate_switchprobs_cM,
    estimate_switchprobs_PS,
)

setup_logging(sm.log[0])

snp_info = sm.input["snp_info"]
tot_mtx_snp = sm.input["tot_mtx_snp"]
a_mtx_snp = sm.input["a_mtx_snp"]
b_mtx_snp = sm.input["b_mtx_snp"]
dp_corrected = sm.input["dp_corrected"]
window_df_file = sm.input["window_df"]

gmap_file = maybe_path(sm.input["gmap_file"])
region_bed = sm.input["region_bed"]
blacklist_bed = maybe_path(sm.input.get("blacklist_bed", None))
genome_size = sm.input["genome_size"]
gtf_file = maybe_path(sm.input["gtf_file"])

qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
run_id = getattr(sm.params, "run_id", "")

sample_df = pd.read_table(sm.input["sample_file"])
sample_name = sample_df["SAMPLE"].iloc[0]
rep_ids = sample_df["REP_ID"].tolist()
sample_types = sample_df["sample_type"].tolist()
assay_type = sm.params["assay_type"]
chromosomes = sm.params["chromosomes"]

has_normal = "normal" in sample_types
tumor_sidx = 1 if has_normal else 0
nsamples = len(sample_df)

logging.info(
    f"combine_counts: sample={sample_name}, assay={assay_type}, "
    f"nsamples={nsamples}, has_normal={has_normal}"
)

snps = pd.read_table(snp_info, sep="\t")
tot_mtx = np.load(tot_mtx_snp)["mat"].astype(np.int32)
a_mtx = np.load(a_mtx_snp)["mat"].astype(np.int32)
b_mtx = np.load(b_mtx_snp)["mat"].astype(np.int32)
nsnps = tot_mtx.shape[0]
logging.info(f"loaded {nsnps} SNPs, {nsamples} samples")

dp_mtx = np.load(dp_corrected)["mat"]
window_df = pd.read_table(window_df_file, sep="\t")
n_window_df = len(window_df)
logging.info(f"loaded {n_window_df} corrected window_df")

assert "region_id" in snps.columns, "invalid SNP file"
grp_cols = ["region_id"]
if "PS" not in snps.columns:
    logging.info("PS not in SNP columns, setting PS=1 for all SNPs")
    snps["PS"] = 1
else:
    logging.info("PS is provided")
num_phaseset = snps["PS"].nunique()
logging.info(f"#phaseset={num_phaseset}")
grp_cols.append("PS")

window_df["win_idx"] = np.arange(len(window_df))
# Assign PS to window_df from SNPs: each window gets the majority-vote PS
_snps_tmp = snps[["#CHR", "POS0", "PS"]].copy()
_snps_tmp = assign_pos_to_range(_snps_tmp, window_df, ref_id="win_idx", pos_col="POS0")
_snps_tmp = _snps_tmp.dropna(subset=["win_idx"])
_snps_tmp["win_idx"] = _snps_tmp["win_idx"].astype(np.int64)
win_ps = _snps_tmp.groupby("win_idx")["PS"].agg(lambda x: x.mode().iloc[0])
window_df["PS"] = window_df["win_idx"].map(win_ps)
# Windows with no SNPs inherit PS from previous window (forward-fill)
if window_df["PS"].isna().any():
    window_df["PS"] = window_df["PS"].ffill()

bbs, snps = adaptive_binning_windows(
    window_df,
    snps,
    tot_mtx,
    int(sm.params["min_snp_reads"]),
    int(sm.params["min_snp_per_block"]),
    grp_cols=grp_cols,
    tumor_sidx=tumor_sidx,
)
num_bbs = len(bbs)
bb_ids = snps["bb_id"].to_numpy()

snp_orig_idx = snps["_orig_idx"].to_numpy()
tot_mtx = tot_mtx[snp_orig_idx]
a_mtx = a_mtx[snp_orig_idx]
b_mtx = b_mtx[snp_orig_idx]

a_mtx_bb = matrix_segmentation(a_mtx, bb_ids, num_bbs)
b_mtx_bb = matrix_segmentation(b_mtx, bb_ids, num_bbs)
tot_mtx_bb = matrix_segmentation(tot_mtx, bb_ids, num_bbs)

plot_allele_freqs(
    bbs,
    rep_ids,
    tot_mtx_bb,
    b_mtx_bb,
    genome_size,
    qc_dir,
    apply_pseudobulk=False,
    allele="B",
    unit="bb",
    region_bed=region_bed,
    blacklist_bed=blacklist_bed,
    run_id=run_id,
)

baf_mtx_bb = np.divide(
    b_mtx_bb,
    tot_mtx_bb,
    where=tot_mtx_bb > 0,
    out=np.full_like(b_mtx_bb, np.nan, dtype=np.float32),
)

logging.info("aggregating corrected window depth into adaptive bins")

win_bin_ids = window_df["bin_id"].to_numpy()
win_lengths = (window_df["END"] - window_df["START"]).to_numpy(dtype=np.float64)
total_len_per_bin = np.bincount(win_bin_ids, weights=win_lengths, minlength=num_bbs)

bb_dp = np.full((num_bbs, nsamples), np.nan, dtype=np.float32)
for s in range(nsamples):
    weighted_sums = np.bincount(
        win_bin_ids, weights=dp_mtx[:, s] * win_lengths, minlength=num_bbs
    )
    with np.errstate(invalid="ignore"):
        bb_dp[:, s] = weighted_sums / total_len_per_bin

log_nan_summary("bb depth", bb_dp, rep_ids, num_bbs)

logging.info(f"compute bb RDR, has_normal={has_normal}")

if has_normal:
    window_sizes = (window_df["END"] - window_df["START"]).to_numpy(dtype=np.float64)
    total_bases = np.nansum(dp_mtx * window_sizes[:, None], axis=0)
    library_correction = total_bases[0] / total_bases[1:]
    logging.info(f"library normalization factor: {library_correction}")

    normal_bb_dp = bb_dp[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        bb_rdr = bb_dp[:, tumor_sidx:] / normal_bb_dp[:, None]
        bb_rdr *= library_correction[None, :]
else:
    bb_rdr = np.full((num_bbs, nsamples), np.nan, dtype=np.float32)
    for i in range(nsamples):
        col = bb_dp[:, i]
        valid_i = np.isfinite(col) & (col > 0)
        if valid_i.any():
            med = np.median(col[valid_i])
            logging.info(f"  bb median-centering {rep_ids[i]}: median={med:.4f}")
            with np.errstate(invalid="ignore", divide="ignore"):
                bb_rdr[valid_i, i] = col[valid_i] / med

n_nan_bb = int(np.isnan(bb_rdr[:, 0]).sum())
_n_filled = num_bbs - n_nan_bb
logging.info(
    f"bb RDR: {_n_filled}/{num_bbs} ({_n_filled / max(num_bbs, 1) * 100:.1f}%) filled, "
    f"{n_nan_bb} NaN"
)

tumor_rep_ids = rep_ids[tumor_sidx:]
bb_rdr_tumor = bb_rdr[:, tumor_sidx:] if not has_normal else bb_rdr
rdr_ylim = (np.round(np.nanquantile(bb_rdr_tumor, 0.99)).astype(int) + 1) * 1.1
plot_1d_multi_sample(
    bbs,
    bb_rdr_tumor,
    list(tumor_rep_ids),
    genome_size,
    stamp_path(os.path.join(qc_dir, "rdr_bb.pdf"), run_id),
    unit="bb",
    val_type="RDR",
    max_ylim=rdr_ylim,
    region_bed=region_bed,
    blacklist_bed=blacklist_bed,
)

nan_mask = (
    np.isnan(baf_mtx_bb).any(axis=1)
    | np.isnan(bb_dp).any(axis=1)
    | np.isnan(bb_rdr).any(axis=1)
)
n_nan_rows = int(nan_mask.sum())
n_valid = num_bbs - n_nan_rows
logging.info(
    f"NaN row filter: {n_nan_rows}/{num_bbs} bins have NaN, "
    f"keeping {n_valid} ({n_valid / max(num_bbs, 1) * 100:.1f}%)"
)

if n_nan_rows > 0:
    valid = ~nan_mask
    bbs = bbs.loc[valid].reset_index(drop=True)
    tot_mtx_bb = tot_mtx_bb[valid]
    a_mtx_bb = a_mtx_bb[valid]
    b_mtx_bb = b_mtx_bb[valid]
    baf_mtx_bb = baf_mtx_bb[valid]
    bb_dp = bb_dp[valid]
    bb_rdr = bb_rdr[valid]

kept_bb_ids = np.where(~nan_mask)[0] if n_nan_rows > 0 else np.arange(num_bbs)
old_to_new = {old: new for new, old in enumerate(kept_bb_ids)}
snps_valid = snps[snps["bb_id"].isin(old_to_new)].copy()
snps_valid["bb_id"] = snps_valid["bb_id"].map(old_to_new)
bbs["bb_id"] = np.arange(len(bbs))

logging.info("estimate bin-level switchprobs")
if gmap_file is not None:
    genetic_map = pd.read_table(gmap_file, sep="\t")
    dist_cms = interp_cM_blocks(bbs, snps_valid, genetic_map, block_id_col="bb_id")
    bbs["switchprobs"] = estimate_switchprobs_cM(
        dist_cms,
        nu=float(sm.params["nu"]),
        min_switchprob=float(sm.params["min_switchprob"]),
    )
else:
    switchprob_ps = float(sm.params["switchprob_ps"])
    bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

bbs.to_csv(sm.output["bb_file"], sep="\t", header=True, index=False)
bbs.to_csv(
    sm.output["bed_file"],
    columns=["#CHR", "START", "END"],
    sep="\t",
    header=False,
    index=False,
)
np.savez_compressed(sm.output["tot_mtx_bb"], mat=tot_mtx_bb)
np.savez_compressed(sm.output["a_mtx_bb"], mat=a_mtx_bb)
np.savez_compressed(sm.output["b_mtx_bb"], mat=b_mtx_bb)
np.savez_compressed(sm.output["baf_mtx_bb"], mat=baf_mtx_bb)
np.savez_compressed(sm.output["dp_mtx_bb"], mat=bb_dp)
np.savez_compressed(sm.output["rdr_mtx_bb"], mat=bb_rdr)
shutil.copy2(sm.input["sample_file"], sm.output["sample_file"])
logging.info("finished combine_counts.")
