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
from scipy import sparse

from utils import read_VCF, symlink_force
from postprocess_utils import *

##################################################
"""
Postprocess bulk data matrices.
"""

logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# inputs
vcf_files = sm.input["vcfs"]
tot_mat_files = sm.input["tot_mats"]
ad_mat_files = sm.input["ad_mats"]
region_bed_file = sm.input["region_bed"]
genome_file = sm.input["genome_size"]

sample_name = sm.params["sample_name"]
rep_ids = sm.params["rep_ids"]

# filtering options
mask_by_region = bool(sm.params["mask_out_of_region"])
min_depth = int(sm.params["min_depth"])
gamma = float(sm.params["gamma"])

##################################################
logging.info(f"postprocess bulk data, sample name={sample_name}")
snps = read_VCF(sm.input["snp_file"], addkey=True)
M = len(snps)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

# add per-SNP subregion info
regions = read_region_file(
    region_bed_file, addchr=str(snps["#CHR"].iloc[0]).startswith("chr")
)
snps = assign_snp_bounderies(snps, regions)
snp_cols = ["#CHR", "POS", "START", "END", "BLOCKSIZE", "region_id"]

# load phase information
gt = snps["GT"].astype(str)
is_phased_arr = gt.str.contains(r"\|", na=False).to_numpy()
is_phased = False
phases = None
if is_phased_arr.any():
    assert is_phased_arr.all(), "some SNPs are unphased"
    snp_cols += ["PS", "GT"]
    is_phased = True
    phases = gt.str[2].astype(np.int8).to_numpy()
snps = snps[snp_cols].reset_index(drop=True)
logging.info(f"SNP file is phased={is_phased}")

# re-order replicate IDs and mat columns.
# ensure rep_ids[0] is "normal" if exists
if "normal" in rep_ids:
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

snp_mask = np.ones(len(snps), dtype=bool)
if mask_by_region:
    snp_mask = snp_mask & get_mask_by_region(snps, regions)
if min_depth > 0:
    snp_mask = snp_mask & get_mask_by_depth(snps, tot_mtx, min_dp=min_depth)
if "normal" in rep_ids:
    snp_mask = snp_mask & get_mask_by_het_balanced(
        snps, ref_mtx, alt_mtx, gamma, normal_idx=0
    )

snps = snps.loc[snp_mask, :].reset_index(drop=True)
tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
alt_mtx = alt_mtx[snp_mask, :]

##################################################
# save outputs
sample_df = pd.DataFrame({"SAMPLE": [f"{sample_name}_{rep_id}" for rep_id in rep_ids]})
sample_df["SAMPLE_NAME"] = sample_name
sample_df["REP_ID"] = rep_ids
sample_df.to_csv(sm.output["sample_file"], sep="\t", header=True, index=False)
snps.to_csv(
    sm.output["bed_file"],
    columns=["#CHR", "START", "END"],
    sep="\t",
    header=False,
    index=False,
)
snps.to_csv(sm.output["info_file"], sep="\t", header=True, index=False)

# save as dense matrix
np.savez_compressed(sm.output["tot_mtx"], mat=tot_mtx)
np.savez_compressed(sm.output["ref_mtx"], mat=ref_mtx)
np.savez_compressed(sm.output["alt_mtx"], mat=alt_mtx)
if is_phased:
    a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, phases[snp_mask])
    a_mtx = a_mtx.astype(np.int32)
    b_mtx = b_mtx.astype(np.int32)
    np.savez_compressed(sm.output["a_mtx"], mat=a_mtx)
    np.savez_compressed(sm.output["b_mtx"], mat=b_mtx)
else:
    symlink_force(sm.output["alt_mtx"], sm.output["a_mtx"])
    symlink_force(sm.output["ref_mtx"], sm.output["b_mtx"])

##################################################
# QC analysis
qc_dir = os.path.join(os.path.dirname(sm.output["tot_mtx"]), "qc")
os.makedirs(qc_dir, exist_ok=True)
plot_snps_allele_freqs(
    snps,
    rep_ids,
    tot_mtx,
    ref_mtx,
    genome_file,
    qc_dir,
    apply_pseudobulk=False,
    allele="ref",
)
if is_phased:
    plot_snps_allele_freqs(
        snps,
        rep_ids,
        tot_mtx,
        b_mtx,
        genome_file,
        qc_dir,
        apply_pseudobulk=False,
        allele="B",
    )

logging.info(f"finished postprocessing bulk data.")
