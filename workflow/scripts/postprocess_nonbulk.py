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

from utils import read_VCF, symlink_force
from postprocess_utils import *

##################################################
"""
Postprocess non-bulk data matrices.
"""

logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# inputs
vcf_files = sm.input["vcfs"]
barcode_files = sm.input["sample_tsvs"]
dp_mat_files = sm.input["dp_mats"]
ad_mat_files = sm.input["ad_mats"]
region_bed_file = sm.input["region_bed"]

sample_name = sm.params["sample_name"]
modality = sm.params["modality"]
rep_ids = list(sm.params["rep_ids"])
data_types = list(sm.params["data_types"])

# filtering options
mask_by_region = bool(sm.params["mask_out_of_region"])

##################################################
logging.info(f"postprocess non-bulk data, sample name={sample_name}")
snps = read_VCF(sm.input["snp_file"], addkey=True)
M = len(snps)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

# load phase information
gt = snps["GT"].astype(str)
is_phased_arr = gt.str.contains(r"\|", na=False).to_numpy()
is_phased = False
phases = None
snp_cols = ["#CHR", "POS"]
if is_phased_arr.any():
    assert is_phased_arr.all(), "some SNPs are unphased"
    snp_cols += ["PS", "GT"]
    is_phased = True
    phases = gt.str[2].astype(np.int8).to_numpy()
snps = snps[snp_cols].reset_index(drop=True)
logging.info(f"SNP file is phased={is_phased}")

snp_mask = np.ones(len(snps), dtype=bool)
if mask_by_region:
    snp_mask = snp_mask & get_mask_by_region(snps, region_bed_file)

if modality == "multiome":
    barcodes_gex = pd.read_table(
        barcode_files[0], sep="\t", header=None, names=["BARCODE"]
    )
    barcodes_atac = pd.read_table(
        barcode_files[1], sep="\t", header=None, names=["BARCODE"]
    )
    assert (barcodes_gex["BARCODE"] == barcodes_atac["BARCODE"]).all()
    barcodes = barcodes_gex
    ncells = len(barcodes)
    rep_id = rep_ids[0]

    dp_mtx_list = []
    ref_mtx_list = []
    alt_mtx_list = []
    for i, data_type in enumerate(data_types):
        dp_mtx, alt_mtx = canon_mat_one_replicate(
            parent_keys, vcf_files[i], dp_mat_files[i], ad_mat_files[i], ncells
        )
        ref_mtx = dp_mtx - alt_mtx
        dp_mtx_list.append(dp_mtx)
        ref_mtx_list.append(ref_mtx)
        alt_mtx_list.append(alt_mtx)
    snps = snps.loc[snp_mask, :].reset_index(drop=True)

    # save outputs
    snps.to_csv(sm.output["info_file"], sep="\t", header=True, index=False)
    snp_ids = snps["#CHR"].astype(str) + "_" + snps["POS"].astype(str)
    np.save(sm.output["unique_snp_ids"], snp_ids.to_numpy())
    barcodes.to_csv(sm.output["barcode_file"], sep="\t", header=False, index=False)
    for i, data_type in enumerate(["scRNA", "scATAC"]):
        dp_mtx = dp_mtx_list[i][snp_mask, :]
        ref_mtx = ref_mtx_list[i][snp_mask, :]
        alt_mtx = alt_mtx_list[i][snp_mask, :]
        # save as sparse matrix
        save_npz(sm.output["dp_mtx"][i], dp_mtx)
        save_npz(sm.output["ref_mtx"][i], ref_mtx)
        save_npz(sm.output["alt_mtx"][i], alt_mtx)
        if is_phased:
            a_mtx, b_mtx = apply_phase_to_mat(
                dp_mtx, ref_mtx, alt_mtx, phases[snp_mask]
            )
            a_mtx = a_mtx.astype(np.int32, copy=False)
            b_mtx = b_mtx.astype(np.int32, copy=False)
            save_npz(sm.output["a_mtx"][i], a_mtx)
            save_npz(sm.output["b_mtx"][i], b_mtx)
        else:
            symlink_force(sm.output["alt_mtx"][i], sm.output["a_mtx"][i])
            symlink_force(sm.output["ref_mtx"][i], sm.output["b_mtx"][i])
else:
    num_samples = len(rep_ids)
    barcode_list = []
    dp_mtx_list = []
    ad_mtx_list = []
    for i, rep_id in enumerate(rep_ids):
        barcodes = pd.read_table(
            barcode_files[i], sep="\t", header=None, names=["BARCODE"]
        )
        if len(rep_ids) > 1:
            barcodes["BARCODE"] = barcodes["BARCODE"].astype(str) + f"_{rep_id}"
        dp_canon, ad_canon = canon_mat_one_replicate(
            parent_keys,
            vcf_files[i],
            dp_mat_files[i],
            ad_mat_files[i],
            ncells=len(barcodes),
        )
        barcode_list.append(barcodes)
        dp_mtx_list.append(dp_canon)
        ad_mtx_list.append(ad_canon)

    all_barcodes = pd.concat(barcode_list, axis=0, ignore_index=True)
    dp_mtx, ref_mtx, alt_mtx = merge_mats(dp_mtx_list, ad_mtx_list)
    snps = snps.loc[snp_mask, :].reset_index(drop=True)

    dp_mtx = dp_mtx[snp_mask, :]
    ref_mtx = ref_mtx[snp_mask, :]
    alt_mtx = alt_mtx[snp_mask, :]

    # save outputs
    snps.to_csv(sm.output["info_file"], sep="\t", header=True, index=False)
    snp_ids = snps["#CHR"].astype(str) + "_" + snps["POS"].astype(str)
    np.save(sm.output["unique_snp_ids"], snp_ids.to_numpy())
    all_barcodes.to_csv(sm.output["all_barcodes"], sep="\t", header=False, index=False)
    
    save_npz(sm.output["dp_mtx"], dp_mtx)
    save_npz(sm.output["ref_mtx"], ref_mtx)
    save_npz(sm.output["alt_mtx"], alt_mtx)
    if is_phased:
        a_mtx, b_mtx = apply_phase_to_mat(dp_mtx, ref_mtx, alt_mtx, phases[snp_mask])
        a_mtx = a_mtx.astype(np.int32, copy=False)
        b_mtx = b_mtx.astype(np.int32, copy=False)
        save_npz(sm.output["a_mtx"], a_mtx)
        save_npz(sm.output["b_mtx"], b_mtx)
    else:
        symlink_force(sm.output["alt_mtx"], sm.output["a_mtx"])
        symlink_force(sm.output["ref_mtx"], sm.output["b_mtx"])

sample_df = pd.DataFrame({"SAMPLE": [f"{sample_name}_{rep_id}" for rep_id in rep_ids]})
sample_df["SAMPLE_NAME"] = sample_name
sample_df["REP_ID"] = rep_ids
sample_df.to_csv(sm.output["sample_file"], sep="\t", header=True, index=False)
logging.info(f"finished postprocessing: #samples={len(rep_ids)}, #SNPs={len(snps)}")
