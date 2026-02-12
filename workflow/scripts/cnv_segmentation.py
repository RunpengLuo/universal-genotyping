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
from scipy.sparse import save_npz, load_npz
import anndata
from scanpy import AnnData
from utils import *
from io_utils import *
from aggregation_utils import *

##################################################
"""
Segment allele/feature-level matrix into CNV segment-level matrix.
Input for Copy-typing.
"""
setup_logging(sm.log[0])

# inputs
snp_info = sm.input["snp_info"]
tot_mtx_snp = sm.input["tot_mtx_snp"]
a_mtx_snp = sm.input["a_mtx_snp"]
b_mtx_snp = sm.input["b_mtx_snp"]
h5ad_file = sm.input["h5ad_file"]

all_barcodes = sm.input["all_barcodes"]
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)

region_bed = sm.input["region_bed"]
genome_size = sm.input["genome_size"]
gtf_file = maybe_path(sm.input["gtf_file"])

sample_name = sm.params["sample_name"]
assay_type = sm.params["assay_type"]
feature_type = sm.params["feature_type"]

sample_df = pd.read_table(sm.input["sample_file"])
rep_ids = sample_df["REP_ID"].tolist()

bulk_assays = {"bulkDNA", "bulkWGS", "bulkWES"}
is_bulk_assay = assay_type in bulk_assays
assert not is_bulk_assay, "bulk sample CNV segmentation unsupported yet"

##################################################
logging.info(f"cnv segmentation, sample name={sample_name}, assay_type={assay_type}")
logging.info(f"rep_ids={rep_ids}")
snps = pd.read_table(snp_info, sep="\t")

if is_bulk_assay:
    # dense matrix
    tot_mtx = np.load(tot_mtx_snp)["mat"].astype(np.int32)
    a_mtx = np.load(a_mtx_snp)["mat"].astype(np.int32)
    b_mtx = np.load(b_mtx_snp)["mat"].astype(np.int32)
else:
    # sparse matrix
    tot_mtx = load_npz(tot_mtx_snp)
    a_mtx = load_npz(a_mtx_snp)
    b_mtx = load_npz(b_mtx_snp)

##################################################
# load cnv profile
segs_df, clones = read_ucn_file(sm.input["seg_ucn"])
num_segs = len(segs_df)
segs_df["seg_id"] = segs_df.index

bbcs_df, _ = read_ucn_file(sm.input["bbc_ucn"])
num_bbcs = len(bbcs_df)
logging.info(f"#CNV segments={num_segs}, #CNV blocks={num_bbcs}")

bbcs_phases_df = pd.read_table(sm.input["bbc_phases"], sep="\t").rename(
    columns={"PHASE": "PHASE-BBC"}
)
bbcs_df = pd.merge(
    left=segs_df,
    right=bbcs_phases_df[["#CHR", "START", "END", "PHASE"]],
    on=["#CHR", "START", "END"],
    how="left",
)
assert bbcs_df["PHASE"].notna().all(), f"corrupted bbc* files, un-matched coordinates"
bbcs_df["bbc_id"] = bbcs_df.index

##################################################
# annotate per-SNP phase by CNV phase
snps["RAW_SNP_IDX"] = np.arange(len(snps))
snps = snp_to_region(snps, bbcs_df, assay_type, region_id="bbc_id")
bbc_phases = pd.merge(
    left=snps, right=bbcs_df[["bbc_id", "PHASE-BBC"]], on="bbc_id", how="left"
)["PHASE-BBC"].to_numpy()

##################################################
# filter mat rows with SNP not present in CNV profile
raw_snp_ids = snps["RAW_SNP_IDX"].to_numpy()
tot_mtx = tot_mtx[raw_snp_ids, :]
a_mtx = a_mtx[raw_snp_ids, :]
b_mtx = b_mtx[raw_snp_ids, :]
a_mtx_corr, b_mtx_corr = apply_phase_to_mat(tot_mtx, b_mtx, a_mtx, bbc_phases)

plot_allele_freqs(
    snps,
    rep_ids,
    tot_mtx,
    b_mtx_corr,
    genome_size,
    qc_dir,
    apply_pseudobulk=not is_bulk_assay,
    allele="cnv-B",
    unit="snp",
    suffix=f"_{assay_type}",
)

# segment counts by CNV profile
snps = snp_to_region(snps, segs_df, assay_type, region_id="seg_id")
seg_ids = snps["seg_id"].to_numpy()
y_count = matrix_segmentation(b_mtx_corr, seg_ids, num_segs)
d_count = matrix_segmentation(tot_mtx, seg_ids, num_segs)
assert y_count.shape[0] == num_segs

##################################################
if not is_bulk_assay:
    adata: AnnData = sc.read_h5ad(h5ad_file)
    # consistent order with allele mats
    barcodes = np.asarray(read_barcodes(all_barcodes), dtype=str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    assert len(missing) == 0, (
        f"Missing {len(missing)} barcodes, e.g. {missing[:5]}, bug!"
    )
    adata = adata[barcodes, :].copy()

    adata = feature_to_blocks(
        adata, segs_df, assay_type, block_idx="seg_id", drop_cols=False
    )
    counts = adata.var["seg_id"].value_counts()
    segs_df[f"#{feature_type}"] = segs_df["seg_id"].map(counts).fillna(0).astype(int)
    x_count = matrix_segmentation(adata.X.T, adata.var["seg_id"].to_numpy(), num_segs)
    save_npz(sm.output["x_count"], x_count)
    save_npz(sm.output["y_count"], y_count)
    save_npz(sm.output["d_count"], d_count)
else:
    np.savez_compressed(sm.output["y_count"], mat=y_count)
    np.savez_compressed(sm.output["d_count"], mat=d_count)

segs_df.to_csv(sm.output["cnv_segments"], header=True, sep="\t", index=False)
logging.info(f"finished.")
