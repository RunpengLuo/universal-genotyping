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
from io_utils import *
from postprocess_utils import *
from aggregate_snps import *
from switchprobs import *

##################################################
logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# inputs
vcf_files = sm.input["vcfs"]
tot_mat_files = sm.input["tot_mats"]
ad_mat_files = sm.input["ad_mats"]
snp_file = sm.input["snp_file"]

region_bed = sm.input["region_bed"]
genome_file = sm.input["genome_size"]
block_bed = maybe_path(sm.input["block_bed"])
gmap_file = sm.input["gmap_file"]

sample_name = sm.params["sample_name"]
rep_ids = sm.params["rep_ids"]

##################################################
logging.info(f"postprocess bulk data, sample name={sample_name}")
snps = read_VCF(snp_file, addkey=True)
snps["POS0"] = snps["POS"] - 1
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

gt = snps["GT"].astype(str)
is_phased_arr = snps["GT"].astype(str).str.contains(r"\|", na=False).to_numpy()
assert is_phased_arr.all(), "some SNPs are unphased"
snps["PHASE"] = gt.str[2].astype(np.int8).to_numpy()

# If normal sample exists, make first column be normal sample data.
has_normal = "normal" in rep_ids
if has_normal:
    normal_idx = list(rep_ids).index("normal")
    rep_ids[0], rep_ids[normal_idx] = rep_ids[normal_idx], rep_ids[0]
    vcf_files[0], vcf_files[normal_idx] = vcf_files[normal_idx], vcf_files[0]
    tot_mat_files[0], tot_mat_files[normal_idx] = (
        tot_mat_files[normal_idx],
        tot_mat_files[0],
    )
    ad_mat_files[0], ad_mat_files[normal_idx] = (
        ad_mat_files[normal_idx],
        ad_mat_files[0],
    )
tumor_sidx = {False: 0, True: 1}[has_normal]


##################################################
# concat to form SNP by sample mats
tot_mtx_list = []
ad_mtx_list = []
for i in range(len(rep_ids)):
    tot_canon, ad_canon = canon_mat_one_replicate(
        parent_keys, vcf_files[i], tot_mat_files[i], ad_mat_files[i], ncells=1
    )
    tot_mtx_list.append(tot_canon)
    ad_mtx_list.append(ad_canon)

tot_mtx, ref_mtx, alt_mtx = merge_mats(tot_mtx_list, ad_mtx_list)
# convert to dense mats for bulk sample
tot_mtx = tot_mtx.toarray()
ref_mtx = ref_mtx.toarray()
alt_mtx = alt_mtx.toarray()

##################################################
# filter invalid SNPs
grp_cols = ["region_id"]
snp_mask = np.ones(len(snps), dtype=bool)
regions = read_region_file(
    region_bed, addchr=str(snps["#CHR"].iloc[0]).startswith("chr")
)
snp_mask = snp_mask & get_mask_by_region(snps, regions)
snp_mask = snp_mask & get_mask_by_depth(
    snps, tot_mtx, min_dp=max(int(sm.params["min_depth"]), 1)
)
if has_normal:
    snp_mask = snp_mask & get_mask_by_het_balanced(
        snps, ref_mtx, alt_mtx, float(sm.params["gamma"]), normal_idx=0
    )

# TODO
# blocks = None
# if block_bed is not None:
#     blocks = read_region_file(
#         block_bed, addchr=str(snps["#CHR"].iloc[0]).startswith("chr")
#     )
#     snps = annotate_snps_to_regions(blocks, snps, seg_id="BLOCK_ID")
#     snp_mask = snp_mask & snps["BLOCK_ID"].notna().to_numpy()
#     grp_cols.append("BLOCK_ID")

snps = snps.loc[snp_mask, :].reset_index(drop=True)
tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
alt_mtx = alt_mtx[snp_mask, :]
(nsnps, nsamples) = tot_mtx.shape

# assign each SNP a consecutive interval [START, END) along chromosome based on midpoints
snps = assign_snp_bounderies(snps, regions, region_id="region_id")

##################################################
# convert to phased A/B mats
a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, snps["PHASE"].to_numpy())
baf_mtx = b_mtx / tot_mtx

qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
plot_allele_freqs(
    snps,
    rep_ids,
    tot_mtx,
    b_mtx,
    genome_file,
    qc_dir,
    apply_pseudobulk=False,
    allele="B",
    unit="snp",
)
# store SNP-level data
snps[["#CHR", "POS", "GT", "PHASE"]].to_csv(
    sm.output["snp_file"], sep="\t", header=True, index=False
)
np.savez_compressed(sm.output["tot_mtx_snp"], mat=tot_mtx)
np.savez_compressed(sm.output["a_mtx_snp"], mat=a_mtx)
np.savez_compressed(sm.output["b_mtx_snp"], mat=b_mtx)


##################################################
# decide aggregation bounderies & estimate switchprobs
if gmap_file is not None:
    logging.info("estimate non-homogeneous switchprobs based on genetic map file")
    genetic_map = pd.read_table(gmap_file, sep="\t")
    dist_cms = interp_cM_snps(snps, genetic_map)
    switchprobs = estimate_switchprobs_cM(
        dist_cms,
        nu=float(sm.params["nu"]),
        min_switchprob=float(sm.params["min_switchprob"]),
    )
    snps["switchprobs"] = switchprobs
    if "PS" not in snps.columns:
        max_switchprob = float(sm.params["max_switchprob"])
        logging.info(f"decide phaseset PS with max_switchprob={max_switchprob}")
        snps["PS"] = (switchprobs > max_switchprob).cumsum()
else:
    assert "PS" in snps.columns, (
        "either PS be provided (long-read) or genetic map should be preovided"
    )
    switchprob_ps = float(sm.params["switchprob_ps"])
    logging.info(f"fix intra-PS switchprobs={switchprob_ps}")
    switchprobs = estimate_switchprobs_PS(snps, switchprob_ps)
    snps["switchprobs"] = switchprobs
grp_cols.append("PS")

binom_test = bool(sm.params["binom_test"])
if binom_test:
    # 2-prop binom test on BAF mats
    switches_2d = np.zeros((nsnps, nsamples), dtype=bool)
    for si in range(tumor_sidx, nsamples):  # tumor samples
        switches_2d[:, si] = binom_2prop_test(
            b_mtx[:, si],
            tot_mtx[:, si],
            baf_mtx[:, si],
            alpha=float(sm.params["binom_alpha"]),
            margin=float(sm.params["binom_margin"]),
        )
    switches = np.any(switches_2d, axis=1)

    # TODO print top 50 BAFs at switched SNPs
    snps["binom_id"] = np.cumsum(switches)
    grp_cols.append("binom_id")

##################################################
# segmentation
bbs = adaptive_binning(
    snps,
    int(sm.params["min_snp_reads"]),
    int(sm.params["min_snp_per_block"]),
    tot_mtx,
    grp_cols,
    colname="bb_id",
    tumor_sidx=tumor_sidx,
)
a_mtx_bb = matrix_segmentation(a_mtx, snps["bb_id"].to_numpy(), len(bbs))
b_mtx_bb = matrix_segmentation(b_mtx, snps["bb_id"].to_numpy(), len(bbs))
tot_mtx_bb = matrix_segmentation(tot_mtx, snps["bb_id"].to_numpy(), len(bbs))
plot_allele_freqs(
    bbs,
    rep_ids,
    tot_mtx_bb,
    b_mtx_bb,
    genome_file,
    qc_dir,
    apply_pseudobulk=False,
    allele="B",
    unit="bb",
)


##################################################
bbs.to_csv(sm.output["bb_file"], sep="\t", header=True, index=False)
np.savez_compressed(sm.output["tot_mtx_bb"], mat=tot_mtx_bb)
np.savez_compressed(sm.output["a_mtx_bb"], mat=a_mtx_bb)
np.savez_compressed(sm.output["b_mtx_bb"], mat=b_mtx_bb)
# RDR position files for mosdepth
bbs.to_csv(
    sm.output["bed_file"],
    columns=["#CHR", "START", "END"],
    sep="\t",
    header=False,
    index=False,
)

sample_df = pd.DataFrame({"SAMPLE": [f"{sample_name}_{rep_id}" for rep_id in rep_ids]})
sample_df["SAMPLE_NAME"] = sample_name
sample_df["REP_ID"] = rep_ids
sample_df.to_csv(sm.output["sample_file"], sep="\t", header=True, index=False)
logging.info(f"finished.")
