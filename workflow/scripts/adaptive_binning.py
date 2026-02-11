import os, sys, gzip, logging
from snakemake.script import snakemake as sm

from collections import OrderedDict

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import save_npz, load_npz

from utils import *
from io_utils import *
from aggregation_utils import *
from postprocess_utils import *

from switchprobs import *

##################################################
"""
Adaptive binning by MSR and MSPB.
Input for HATCHet3 and CalicoST.
"""
setup_logging(sm.log[0])

snp_file = sm.input["snp_file"]
tot_mtx_snp = sm.input["tot_mtx_snp"]
a_mtx_snp = sm.input["a_mtx_snp"]
b_mtx_snp = sm.input["b_mtx_snp"]

gmap_file = maybe_path(sm.input["gmap_file"])
all_barcodes = sm.input["all_barcodes"]
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)

region_bed = sm.input["region_bed"]
genome_size = sm.input["genome_size"]
gtf_file = maybe_path(sm.input["gtf_file"])

sample_df = pd.read_table(sm.input["sample_file"])
sample_name = sample_df["SAMPLE"].iloc[0]
rep_ids = sample_df["REP_ID"].tolist()
assay_type = sm.params["assay_type"]

bulk_assays = {"bulkDNA", "bulkWGS", "bulkWES"}
is_bulk_assay = assay_type in bulk_assays
has_normal = "normal" in rep_ids
tumor_sidx = {False: 0, True: 1}[has_normal]

##################################################
logging.info(
    f"adaptive binning, sample name={sample_name}, assay_type={assay_type}, is_bulk_assay={is_bulk_assay}"
)
logging.info(f"rep_ids={rep_ids}")
snps = pd.read_table(snp_file, sep="\t")
# barcodes = pd.read_table(all_barcodes, sep="\t", header=None, names=["BARCODE"])

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
(nsnps, nsamples) = tot_mtx.shape

##################################################
# decide aggregation bounderies & estimate switchprobs
assert "region_id" in snps.columns, "invalid SNP file"
grp_cols = ["region_id"]

# 1. estimate phaseset if N/A
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
num_phaseset = snps["PS"].nunique()
logging.info(f"#phaseset={num_phaseset}")
grp_cols.append("PS")

# 2. estimate phase error bounderies
binom_test = bool(sm.params["binom_test"])
if binom_test:
    if is_bulk_assay:
        baf_mtx = b_mtx / tot_mtx
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

        num_switches = np.sum(switches)
        logging.info(f"#binom_test switches={num_switches}/{len(switches) - 1}")
        if num_switches > 0:
            pair_mask = switches[1:]  # N-1
            prev_bafs = baf_mtx[:, tumor_sidx:][:-1, :][pair_mask, :]
            next_bafs = baf_mtx[:, tumor_sidx:][1:, :][pair_mask, :]
            avg_abs_bafdevs = np.mean(np.abs(next_bafs - prev_bafs), axis=1)
            edges = np.linspace(0.0, 1.0, 21)
            counts, _ = np.histogram(avg_abs_bafdevs, bins=edges)
            logging.info("pairwise SNP average BAF absolute deviations")
            logging.info("bin_left\tbin_right\tcount")
            for l, r, c in zip(edges[:-1], edges[1:], counts):
                logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")
            snps["binom_id"] = np.cumsum(switches)
            grp_cols.append("binom_id")
    else:
        # binomial test on pseudobulk sample,
        # TODO use spatial information to use tumor-like spots?
        b_vec = np.asarray(b_mtx.sum(axis=1)).ravel()
        tot_vec = np.asarray(tot_mtx.sum(axis=1)).ravel()
        baf_vec = compute_af_pseudobulk(tot_mtx, b_mtx)
        switches = binom_2prop_test(
            b_vec,
            tot_vec,
            baf_vec,
            alpha=float(sm.params["binom_alpha"]),
            margin=float(sm.params["binom_margin"]),
        )
        num_switches = np.sum(switches)
        logging.info(f"#binom_test switches={num_switches}/{len(switches) - 1}")
        if num_switches > 0:
            pair_mask = switches[1:]  # N-1
            prev_bafs = baf_vec[:-1][pair_mask]
            next_bafs = baf_vec[1:][pair_mask]
            avg_abs_bafdevs = np.abs(next_bafs - prev_bafs)
            edges = np.linspace(0.0, 1.0, 21)
            counts, _ = np.histogram(avg_abs_bafdevs, bins=edges)
            logging.info("pairwise SNP average BAF absolute deviations")
            logging.info("bin_left\tbin_right\tcount")
            for l, r, c in zip(edges[:-1], edges[1:], counts):
                logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")
            snps["binom_id"] = np.cumsum(switches)
            grp_cols.append("binom_id")

# 3. use HMM to estimate bounderies? TODO
# 4. restrict within gene binning by feature_id? TODO
##################################################
# meta-snp segmentation
meta_snps = adaptive_binning(
    snps,
    0,
    int(sm.params["nsnp_meta"]),
    tot_vec[:, None],
    grp_cols,
    colname="meta_id",
)
meta_ids = snps["meta_id"].to_numpy()
a_mtx_meta = matrix_segmentation(a_mtx, meta_ids, len(meta_snps))
b_mtx_meta = matrix_segmentation(b_mtx, meta_ids, len(meta_snps))
tot_mtx_meta = matrix_segmentation(tot_mtx, meta_ids, len(meta_snps))

plot_allele_freqs(
    meta_snps,
    rep_ids,
    tot_mtx_meta,
    b_mtx_meta,
    genome_size,
    qc_dir,
    apply_pseudobulk=not is_bulk_assay,
    allele="B",
    unit="meta-snp",
)
meta_snps.to_csv(sm.params["meta_file"], sep="\t", header=True, index=False)
save_npz(sm.params["tot_mtx_meta"], tot_mtx_meta)
save_npz(sm.params["a_mtx_meta"], a_mtx_meta)
save_npz(sm.params["b_mtx_meta"], b_mtx_meta)

##################################################
# MSR&MSPB block segmentation
bbs = adaptive_binning(
    snps,
    int(sm.params["min_snp_reads"]),
    int(sm.params["min_snp_per_block"]),
    tot_mtx,
    grp_cols,
    colname="bb_id",
    tumor_sidx=tumor_sidx,
)

bb_ids = snps["bb_id"].to_numpy()
num_bbs = len(bbs)
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
    apply_pseudobulk=not is_bulk_assay,
    allele="B",
    unit="bb",
)

bbs.to_csv(sm.output["bb_file"], sep="\t", header=True, index=False)
# RDR position files for mosdepth
bbs.to_csv(
    sm.output["bed_file"],
    columns=["#CHR", "START", "END"],
    sep="\t",
    header=False,
    index=False,
)
if is_bulk_assay:
    np.savez_compressed(sm.output["tot_mtx_bb"], mat=tot_mtx_bb)
    np.savez_compressed(sm.output["a_mtx_bb"], mat=a_mtx_bb)
    np.savez_compressed(sm.output["b_mtx_bb"], mat=b_mtx_bb)
else:
    save_npz(sm.params["tot_mtx_bb"], tot_mtx_bb)
    save_npz(sm.params["a_mtx_bb"], a_mtx_bb)
    save_npz(sm.params["b_mtx_bb"], b_mtx_bb)
logging.info(f"finished.")
