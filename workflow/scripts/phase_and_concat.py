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
from scipy.sparse import save_npz

from utils import *
from io_utils import *
from postprocess_utils import *
from aggregation_utils import *
from switchprobs import *

##################################################
"""
Phase and concat allele-level count matrices.
"""
setup_logging(sm.log[0])

# inputs
vcf_files = sm.input["vcfs"]
sample_tsvs = sm.input["sample_tsvs"]
tot_mtx_files = sm.input["tot_mtxs"]
ad_mtx_files = sm.input["ad_mtxs"]
snp_file = sm.input["snp_file"]
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)

region_bed = sm.input["region_bed"]
genome_size = sm.input["genome_size"]
gtf_file = maybe_path(sm.input["gtf_file"])

sample_name = sm.params["sample_name"]
assay_type = sm.params["assay_type"]
rep_ids = list(sm.params["rep_ids"])

bulk_assays = {"bulkWGS", "bulkWES"}
is_bulk_assay = assay_type in bulk_assays

##################################################
logging.info(
    f"phase and concat, sample name={sample_name}, assay_type={assay_type}, is_bulk_assay={is_bulk_assay}"
)
logging.info(f"rep_ids={rep_ids}")

snps = read_VCF(snp_file, addkey=True, add_phase1=True, add_pos0=True)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

if is_bulk_assay:
    has_normal = "normal" in rep_ids
    # If normal sample exists, make first column be normal sample.
    if has_normal:
        normal_idx = list(rep_ids).index("normal")
        rep_ids[0], rep_ids[normal_idx] = rep_ids[normal_idx], rep_ids[0]
        vcf_files[0], vcf_files[normal_idx] = vcf_files[normal_idx], vcf_files[0]
        tot_mtx_files[0], tot_mtx_files[normal_idx] = (
            tot_mtx_files[normal_idx],
            tot_mtx_files[0],
        )
        ad_mtx_files[0], ad_mtx_files[normal_idx] = (
            ad_mtx_files[normal_idx],
            ad_mtx_files[0],
        )

# concat to form SNP by sample mats
barcodes_list = []
tot_mtx_list = []
ad_mtx_list = []
for idx, rep_id in enumerate(rep_ids):
    barcodes = pd.read_table(sample_tsvs[idx], sep="\t", header=None, names=["BARCODE"])
    barcodes["BARCODE"] = barcodes["BARCODE"].astype(str) + f"_{rep_id}"
    nbarcodes = len(barcodes)
    barcodes_list.append(barcodes)
    if is_bulk_assay:
        nbarcodes = 1
    tot_canon, ad_canon = canon_mat_one_replicate(
        parent_keys, vcf_files[idx], tot_mtx_files[idx], ad_mtx_files[idx], nbarcodes
    )
    tot_mtx_list.append(tot_canon)
    ad_mtx_list.append(ad_canon)

all_barcodes = pd.concat(barcodes_list, axis=0, ignore_index=True)
tot_mtx, ref_mtx, alt_mtx = merge_mats(tot_mtx_list, ad_mtx_list)
a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, snps["PHASE"].to_numpy())

# convert to dense mats for bulk sample
if is_bulk_assay:
    tot_mtx = tot_mtx.toarray()
    ref_mtx = ref_mtx.toarray()
    alt_mtx = alt_mtx.toarray()
    a_mtx = a_mtx.toarray()
    b_mtx = b_mtx.toarray()

##################################################
# annotate&filter SNPs
num_snps_before = len(snps)
# define continguous autosome segments
regions = read_region_file(region_bed)
if is_bulk_assay:
    snp_mask = np.ones(len(snps), dtype=bool)
    snp_mask = snp_mask & get_mask_by_region(snps, regions)
    snp_mask = snp_mask & get_mask_by_depth(
        snps, tot_mtx, min_dp=max(int(sm.params["min_depth"]), 1)
    )
    if has_normal:
        snp_mask = snp_mask & get_mask_by_het_balanced(
            snps, ref_mtx, alt_mtx, float(sm.params["gamma"]), normal_idx=0
        )

    if assay_type == "bulkWES":
        logging.info(
            "TODO, filter non-exonic SNPs, restrict SNPs interval by exonic region."
        )
        pass
    snps = snps.loc[snp_mask, :].reset_index(drop=True)
    snps["feature_id"] = "Unknown"
    snps = assign_snp_bounderies(snps, regions, colname="region_id")
else:
    logging.info(f"annotate SNPs with feature_id")
    adata: sc.AnnData = sc.read_h5ad(sm.input["h5ad_file"])
    feature_df = adata.var.reset_index(drop=False).rename(
        columns={"index": "feature_id"}
    )
    feature_df["feature_idx"] = np.arange(len(feature_df))

    snps = assign_pos_to_range(snps, feature_df, ref_id="feature_idx", pos_col="POS0")
    snp_mask = snps["feature_idx"].notna().to_numpy()
    logging.info(
        f"#SNPs overlap with {assay_type} features={np.sum(snp_mask) / len(snps):.3%}"
    )
    snps = snps.loc[snp_mask, :].reset_index(drop=True)
    snps["feature_idx"] = snps["feature_idx"].astype(feature_df["feature_idx"].dtype)
    snps = pd.merge(
        left=snps,
        right=feature_df[["feature_idx", "feature_id"]],
        how="left",
        on="feature_idx",
        sort=False,
    ).reset_index(drop=True)

    snps = assign_snp_bounderies(snps, regions, colname="region_id")
    # fake START/END field, TODO refine by gene interval?
    snps["START"] = snps["POS0"]
    snps["END"] = snps["POS"]

num_snps_after = np.sum(snp_mask)
logging.info(f"#SNPs={num_snps_after}/{num_snps_before} after SNP filtering")

tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
alt_mtx = alt_mtx[snp_mask, :]
a_mtx = a_mtx[snp_mask, :]
b_mtx = b_mtx[snp_mask, :]


##################################################
# QC plots
plot_allele_freqs(
    snps,
    rep_ids,
    tot_mtx,
    ref_mtx,
    genome_size,
    qc_dir,
    apply_pseudobulk=not is_bulk_assay,
    allele="ref",
    unit="snp",
    suffix=f"_{assay_type}",
)
plot_allele_freqs(
    snps,
    rep_ids,
    tot_mtx,
    b_mtx,
    genome_size,
    qc_dir,
    apply_pseudobulk=not is_bulk_assay,
    allele="B",
    unit="snp",
    suffix=f"_{assay_type}",
)

##################################################
logging.info(f"saving output files")
snps[["#CHR", "POS", "START", "END", "GT", "PHASE", "region_id", "feature_id"]].to_csv(
    sm.output["snp_file"], sep="\t", header=True, index=False
)
if is_bulk_assay:
    np.savez_compressed(sm.output["tot_mtx_snp"], mat=tot_mtx)
    np.savez_compressed(sm.output["a_mtx_snp"], mat=a_mtx)
    np.savez_compressed(sm.output["b_mtx_snp"], mat=b_mtx)
else:
    save_npz(sm.output["tot_mtx_snp"], tot_mtx)
    save_npz(sm.output["a_mtx_snp"], a_mtx)
    save_npz(sm.output["b_mtx_snp"], b_mtx)
    # legacy CalicoST allele-level data
    snp_ids = snps["#CHR"].astype(str) + "_" + snps["POS"].astype(str)
    np.save(sm.output["unique_snp_ids"], snp_ids.to_numpy())
    save_npz(sm.output["cell_snp_Aallele"], a_mtx)
    save_npz(sm.output["cell_snp_Ballele"], b_mtx)
all_barcodes.to_csv(sm.output["all_barcodes"], sep="\t", header=False, index=False)
sample_df = pd.DataFrame({"SAMPLE": [f"{sample_name}_{rep_id}" for rep_id in rep_ids]})
sample_df["SAMPLE_NAME"] = sample_name
sample_df["REP_ID"] = rep_ids
sample_df.to_csv(sm.output["sample_file"], sep="\t", header=True, index=False)
logging.info(f"finished.")
