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
from switchprobs import *
from aggregate_snps import *

##################################################
"""
Handle one replicate of multiome data or multiple replicates of single-modality data.
"""
logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# inputs
vcf_files = sm.input["vcfs"]
barcode_files = sm.input["sample_tsvs"]
tot_mat_files = sm.input["tot_mats"]
ad_mat_files = sm.input["ad_mats"]
snp_file = sm.input["snp_file"]

region_bed = sm.input["region_bed"]
genome_size = sm.input["genome_size"]
block_bed = maybe_path(sm.input["block_bed"])

sample_name = sm.params["sample_name"]
modality = sm.params["modality"]
rep_ids = list(sm.params["rep_ids"])
data_types = list(sm.params["data_types"])
has_cn_profile = bool(sm.params["has_cn_profile"])

##################################################
logging.info(
    f"postprocess {modality} data, sample name={sample_name}, data_types={data_types}, rep_ids={rep_ids}"
)
snps = read_VCF(snp_file, addkey=True)
snps["POS0"] = snps["POS"] - 1
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

gt = snps["GT"].astype(str)
is_phased_arr = snps["GT"].astype(str).str.contains(r"\|", na=False).to_numpy()
assert is_phased_arr.all(), "some SNPs are unphased"
snps["PHASE"] = gt.str[2].astype(np.int8).to_numpy()
phases = snps["PHASE"].to_numpy()

##################################################
# concat to form SNP by barcode mats
tot_mtxs = []
ref_mtxs = []
alt_mtxs = []
if modality == "multiome":
    barcodes_gex = pd.read_table(
        barcode_files[0], sep="\t", header=None, names=["BARCODE"]
    )
    barcodes_atac = pd.read_table(
        barcode_files[1], sep="\t", header=None, names=["BARCODE"]
    )
    assert (barcodes_gex["BARCODE"] == barcodes_atac["BARCODE"]).all()
    all_barcodes = barcodes_gex
    nbarcodes = len(all_barcodes)
    for i, data_type in enumerate(data_types):
        tot_mtx, alt_mtx = canon_mat_one_replicate(
            parent_keys, vcf_files[i], tot_mat_files[i], ad_mat_files[i], nbarcodes
        )
        ref_mtx = tot_mtx - alt_mtx
        tot_mtxs.append(tot_mtx)
        ref_mtxs.append(ref_mtx)
        alt_mtxs.append(alt_mtx)
else:
    barcode_list = []
    tot_mtx_list = []
    ad_mtx_list = []
    for i, rep_id in enumerate(rep_ids):
        barcodes = pd.read_table(
            barcode_files[i], sep="\t", header=None, names=["BARCODE"]
        )
        if len(rep_ids) > 1:
            barcodes["BARCODE"] = barcodes["BARCODE"].astype(str) + f"_{rep_id}"
        tot_canon, ad_canon = canon_mat_one_replicate(
            parent_keys,
            vcf_files[i],
            tot_mat_files[i],
            ad_mat_files[i],
            ncells=len(barcodes),
        )
        barcode_list.append(barcodes)
        tot_mtx_list.append(tot_canon)
        ad_mtx_list.append(ad_canon)

    all_barcodes = pd.concat(barcode_list, axis=0, ignore_index=True)
    tot_mtx, ref_mtx, alt_mtx = merge_mats(tot_mtx_list, ad_mtx_list)
    tot_mtxs = [tot_mtx]
    ref_mtxs = [ref_mtx]
    alt_mtxs = [alt_mtx]

##################################################
# apply SNP-level phases
a_mtxs = []
b_mtxs = []
for i in range(len(data_types)):
    a_mtx, b_mtx = apply_phase_to_mat(tot_mtxs[i], ref_mtxs[i], alt_mtxs[i], phases)
    a_mtxs.append(a_mtx)
    b_mtxs.append(b_mtx)

if modality in ["VISIUM", "VISIUM3prime"]:
    logging.info("store legacy CalicoST SNP-level phased mats for VISIUM* data")
    snp_ids = snps["#CHR"].astype(str) + "_" + snps["POS"].astype(str)
    np.save(sm.output["unique_snp_ids_legacy"], snp_ids.to_numpy())
    save_npz(sm.output["a_mtx_legacy"], a_mtxs[0])
    save_npz(sm.output["b_mtx_legacy"], b_mtxs[0])

##################################################
# filter invalid SNPs
grp_cols = ["region_id"]
snp_mask = np.ones(len(snps), dtype=bool)
regions = read_region_file(
    region_bed, addchr=str(snps["#CHR"].iloc[0]).startswith("chr")
)
snp_mask = snp_mask & get_mask_by_region(snps, regions)

# TODO how to handle gene blocks? check numbat-multiome, calicoST has 0 doc about this.
# blocks = None
# if block_bed is not None:
#     blocks = read_region_file(
#         block_bed, addchr=str(snps["#CHR"].iloc[0]).startswith("chr")
#     )
#     snps = annotate_snps_to_regions(blocks, snps, seg_id="BLOCK_ID")
#     snp_mask = snp_mask & snps["BLOCK_ID"].notna().to_numpy()
#     grp_cols.append("BLOCK_ID")

snps = snps.loc[snp_mask, :].reset_index(drop=True)
for i in range(len(data_types)):
    tot_mtxs[i] = tot_mtxs[i][snp_mask, :]
    ref_mtxs[i] = ref_mtxs[i][snp_mask, :]
    alt_mtxs[i] = alt_mtxs[i][snp_mask, :]
    a_mtxs[i] = a_mtxs[i][snp_mask, :]
    b_mtxs[i] = b_mtxs[i][snp_mask, :]

##################################################
# QC - SNP-level phases
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
for i, data_type in enumerate(data_types):
    a_mtx, b_mtx, tot_mtx = a_mtxs[i], b_mtxs[i], tot_mtxs[i]
    plot_allele_freqs(
        snps,
        rep_ids,
        tot_mtx,
        b_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        allele="B",
        unit="snp",
        suffix=f"_{data_type}",
    )
    np.savez_compressed(sm.output["tot_mtx_snp"][i], mat=tot_mtx)
    np.savez_compressed(sm.output["a_mtx_snp"][i], mat=a_mtx)
    np.savez_compressed(sm.output["b_mtx_snp"][i], mat=b_mtx)

# store SNP-level data
snps[["#CHR", "POS", "GT", "PHASE"]].to_csv(
    sm.output["snp_file"], sep="\t", header=True, index=False
)

if has_cn_profile:
    ##################################################
    segs_df, clones, clone_props = read_ucn_file(sm.input["seg_ucn"])
    num_segs = len(segs_df)
    bbcs_df, _, _ = read_ucn_file(sm.input["bbc_ucn"])
    bbcs_phases_df = pd.read_table(sm.input["bbc_phases"], sep="\t").rename(
        columns={"PHASE": "PHASE-BBC"}
    )
    bbcs_df = pd.merge(
        left=segs_df,
        right=bbcs_phases_df[["#CHR", "START", "END", "PHASE"]],
        on=["#CHR", "START", "END"],
        how="left",
    )
    assert bbcs_df["PHASE"].notna().all(), (
        f"corrupted bbc* files, un-matched coordinates"
    )
    bbcs_df["bbc_id"] = bbcs_df.index

    for i, data_type in enumerate(data_types):
        tot_mtx, ref_mtx, alt_mtx = tot_mtxs[i], ref_mtxs[i], alt_mtxs[i]
        (
            dt_snps,
            tot_mtx,
            ref_mtx,
            alt_mtx,
            a_mtx,
            b_mtx,
            segs_df,
            _,
            y_count,
            d_count,
        ) = cnv_guided_allele_aggregation(
            snps, segs_df, bbcs_df, data_type, tot_mtx, ref_mtx, alt_mtx, "bbc_id"
        )
        plot_allele_freqs(
            dt_snps,
            rep_ids,
            tot_mtx,
            b_mtx,
            genome_size,
            qc_dir,
            apply_pseudobulk=True,
            allele="cnv-B",
            unit="snp",
            suffix=f"_{data_type}",
        )
        save_npz(sm.output["y_count"][i], y_count)
        save_npz(sm.output["d_count"][i], d_count)
    segs_df.to_csv(sm.output["cnp_file"], header=True, sep="\t", index=False)
else:
    ##################################################
    # TODO, for multiome data, block_bed file to decide common blocks for gene/peak?
    assert modality != "multiome", "adaptive binning for multiome data TODO"
    # TODO assign intervals by block_bed? e.g., evenly divides gencode ranges.
    snps = assign_snp_bounderies(snps, regions, region_id="region_id")

    a_mtx, b_mtx, tot_mtx = a_mtxs[0], b_mtxs[0], tot_mtxs[0]
    # pseudobulk count vectors
    b_vec = np.asarray(b_mtx.sum(axis=1)).ravel()
    tot_vec = np.asarray(tot_mtx.sum(axis=1)).ravel()
    ##################################################
    # decide aggregation bounderies
    genetic_map = pd.read_table(sm.input["gmap_file"], sep="\t")
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
    grp_cols.append("PS")

    binom_test = bool(sm.params["binom_test"])
    if binom_test:
        # binomial test on pseudobulk sample, TODO use spatial information to use tumor-like spots?
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
        apply_pseudobulk=True,
        allele="B",
        unit="meta-snp",
    )
    meta_snps.to_csv(sm.output["meta_file"][0], sep="\t", header=True, index=False)
    save_npz(sm.output["tot_mtx_meta"][0], tot_mtx_meta)
    save_npz(sm.output["a_mtx_meta"][0], a_mtx_meta)
    save_npz(sm.output["b_mtx_meta"][0], b_mtx_meta)

    # MSR&MSPB segmentation
    bbs = adaptive_binning(
        snps,
        int(sm.params["min_snp_reads"]),
        int(sm.params["min_snp_per_block"]),
        tot_vec[:, None],
        grp_cols,
        colname="bb_id",
    )
    bb_ids = snps["bb_id"].to_numpy()
    a_mtx_bb = matrix_segmentation(a_mtx, bb_ids, len(bbs))
    b_mtx_bb = matrix_segmentation(b_mtx, bb_ids, len(bbs))
    tot_mtx_bb = matrix_segmentation(tot_mtx, bb_ids, len(bbs))

    plot_allele_freqs(
        bbs,
        rep_ids,
        tot_mtx_bb,
        b_mtx_bb,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        allele="B",
        unit="haplo-block",
    )

    bbs.to_csv(sm.output["bb_file"][0], sep="\t", header=True, index=False)
    save_npz(sm.output["tot_mtx_bb"][0], tot_mtx_bb)
    save_npz(sm.output["a_mtx_bb"][0], a_mtx_bb)
    save_npz(sm.output["b_mtx_bb"][0], b_mtx_bb)

sample_df = pd.DataFrame({"SAMPLE": [f"{sample_name}_{rep_id}" for rep_id in rep_ids]})
sample_df["SAMPLE_NAME"] = sample_name
sample_df["REP_ID"] = rep_ids
sample_df.to_csv(sm.output["sample_file"], sep="\t", header=True, index=False)
all_barcodes.to_csv(sm.output["all_barcodes"], sep="\t", header=False, index=False)
logging.info(f"finished.")
