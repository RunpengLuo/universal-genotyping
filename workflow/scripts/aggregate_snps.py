import os
import logging

import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy.sparse import csr_matrix, hstack, issparse
from scipy.stats import beta

import pyranges as pr

from io_utils import *
from postprocess_utils import *


def assign_pos_to_range(
    qry: pd.DataFrame,
    ref: pd.DataFrame,
    ref_id="SUPER_VAR_IDX",
    pos_col="POS0",
    nodup=True,
):
    qry_pr = pr.PyRanges(
        chromosomes=qry["#CHR"], starts=qry[pos_col], ends=qry[pos_col] + 1
    )
    qry_pr = qry_pr.insert(pd.Series(data=qry.index.to_numpy(), name="qry_index"))
    ref_pr = pr.PyRanges(
        chromosomes=ref["#CHR"],
        starts=ref["START"],
        ends=ref["END"],
    )
    ref_pr = ref_pr.insert(ref[ref_id])

    if nodup:  # qry assigned to at most one ref.
        joined = qry_pr.join(ref_pr)
        qry[ref_id] = pd.NA
        qry.loc[joined.df["qry_index"], ref_id] = joined.df[ref_id].to_numpy()
        return qry
    else:
        hits = (
            qry_pr.join(ref_pr)
            .df[["Chromosome", "Start", "qry_index", ref_id]]
            .rename(columns={"Chromosome": "#CHR", "Start": "POS0"})
        )
        hits["POS"] = hits["POS0"] + 1
        return hits


def snp_to_region(
    snp_df: pd.DataFrame, region_df: pd.DataFrame, modality: str, region_id="BIN_ID"
):
    """
    region_df must be 0-indexed non-overlapping intervals [s, t) in standard BED format.
    """
    logging.info(f"#{modality}-SNP (raw)={len(snp_df)}")
    snp_df = assign_pos_to_range(snp_df, region_df, ref_id=region_id, pos_col="POS0")
    isna_snp_df = snp_df[region_id].isna()
    logging.info(
        f"#{modality}-SNPs outside any region={np.sum(isna_snp_df) / len(snp_df):.3%}"
    )
    snp_df.dropna(subset=region_id, inplace=True)
    snp_df[region_id] = snp_df[region_id].astype(region_df[region_id].dtype)
    logging.info(f"#{modality}-SNP (remain)={len(snp_df)}")

    counts = snp_df[region_id].value_counts()
    region_df["#SNPS"] = region_df[region_id].map(counts).fillna(0).astype(int)
    return snp_df


def adaptive_binning(
    snps: pd.DataFrame,
    min_snp_reads: int,
    min_snp_per_block: int,
    tot_mtx: np.ndarray,
    grp_cols: list,
    colname="",
    tumor_sidx=0,
):
    """
    Adaptive binning over pre-defined chromosome regions,
    each bin has MSR and MSPB constraint, beside last bin in the region.
    in-place add <colname> column to snps to indicate bin ids.
    """
    bin_id = 0
    snps[colname] = 0
    snp_grps = snps.groupby(by=grp_cols, sort=False)
    logging.info(f"#SNP groups={len(snp_grps)}, grouper: {grp_cols}")
    for _, grp_snps in snp_grps:
        grp_idxs = grp_snps.index.to_numpy()
        nsnp_grp = len(grp_idxs)
        bin_ids = np.zeros(nsnp_grp, dtype=np.int64)
        bin_id0 = bin_id
        prev_start = 0
        acc_read_count = tot_mtx[grp_idxs[0], tumor_sidx:].copy()
        acc_num_snp = 1
        for i in range(1, nsnp_grp):
            idx = grp_idxs[i]
            if (
                np.all(acc_read_count >= min_snp_reads)
                and acc_num_snp >= min_snp_per_block
            ):
                bin_ids[prev_start:i] = bin_id
                bin_id += 1
                prev_start = i

                acc_read_count = tot_mtx[idx, tumor_sidx:].copy()
                acc_num_snp = 1
            else:
                # extend feature block
                acc_read_count += tot_mtx[idx, tumor_sidx:].copy()
                acc_num_snp += 1

        # fill last block if any
        bin_ids[prev_start:] = max(bin_id - 1, bin_id0)
        snps.loc[grp_idxs, bin_id] = bin_ids

        # if only a partial block is found, block_id is not incremented in the loop
        bin_id = max(bin_id, bin_id0 + 1)

    pos_dict = {
        "#CHR": ("#CHR", "first"),
        "START": ("START", "min"),
        "END": ("END", "max"),
        # effective first/last SNP positions.
        "START0": ("POS0", "min"),
        "END0": ("POS", "max"),
        # inter-block switchprobs
        "switchprobs": ("switchprobs", "first"),
    }
    for grp_col in grp_cols:
        pos_dict[grp_col] = (grp_col, "first")

    snp_grps = snps.groupby(by=bin_id, sort=False, as_index=True)
    snp_bins = snp_grps.agg(**pos_dict)
    snp_bins.loc[:, "#SNPS"] = snp_grps.size().reset_index(drop=True)
    snp_bins.loc[:, "BLOCKSIZE"] = snp_bins["END"] - snp_bins["START"]
    snp_bins[bin_id] = snp_bins.index
    return snp_bins


def matrix_segmentation(X, bin_ids, K):
    """

    N: #snps
    M: #samples
    X: (N, M) sparse or dense   [snp-by-sample]
    bin_ids: (N,) ints in [0..K-1]  (assign each SNP to a bin, bin ids are non-decreasing)
    return:
      - sparse in  -> (K, M) csr_matrix
      - dense in   -> (K, M) ndarray
    """
    sparse_in = issparse(X)
    X = X.tocsr() if sparse_in else np.asarray(X)

    bin_ids = np.asarray(bin_ids, dtype=np.int64)
    N, M = X.shape
    if bin_ids.shape[0] != N:
        raise ValueError(f"bin_ids length {bin_ids.shape[0]} != N {N}")
    if N and (bin_ids.min() < 0 or bin_ids.max() >= K):
        raise ValueError("bin_ids out of range")

    # (K, N) one-hot: row=bin, col=snp
    B = csr_matrix(
        (np.ones(N, dtype=np.int8), (bin_ids, np.arange(N, dtype=np.int64))),
        shape=(K, N),
    )

    X_bin = B @ X  # (K, M)

    if sparse_in:
        return X_bin.tocsr()
    else:
        return X_bin.toarray()  # ensure dense output


def cnv_guided_allele_aggregation(
    snps: pd.DataFrame,
    segs_df: pd.DataFrame,
    bbcs_df: pd.DataFrame,
    data_type: str,
    tot_mtx=None,
    ref_mtx=None,
    alt_mtx=None,
    region_id="bbc_id",
):
    """
    1. integrate HATCHet phases to correct initial phase.
    2. build phased mats
    """
    snps = snps.copy(deep=True)
    # map SNPs to CNV bins, align with HMM clustering step of HATCHet
    snps["RAW_SNP_IDX"] = np.arange(len(snps))
    snps = snp_to_region(snps, bbcs_df, data_type, region_id=region_id)
    phase_inds = pd.merge(
        left=snps, right=bbcs_df[[region_id, "PHASE-BBC"]], on=region_id, how="left"
    )["PHASE-BBC"].to_numpy()
    raw_phase = snps["PHASE"].to_numpy()
    corr_phase = raw_phase * phase_inds + (1 - raw_phase) * (1 - phase_inds)
    snps["PHASE-CORR"] = corr_phase

    # filter mat rows with SNP not present in HATCHet CNV profile
    raw_snp_ids = snps["RAW_SNP_IDX"].to_numpy()
    tot_mat = tot_mat[raw_snp_ids, :]
    ref_mat = ref_mat[raw_snp_ids, :]
    alt_mat = alt_mat[raw_snp_ids, :]
    a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, corr_phase)

    # aggregate SNP-level mats to CNV segment-level mats
    snps = snp_to_region(snps, segs_df, data_type, region_id="seg_id")
    seg_ids = snps["seg_id"].to_numpy()
    num_segs = len(segs_df)
    a_mtx_seg = matrix_segmentation(a_mtx, seg_ids, num_segs)
    b_mtx_seg = matrix_segmentation(b_mtx, seg_ids, num_segs)
    tot_mtx_seg = matrix_segmentation(tot_mtx, seg_ids, num_segs)
    assert a_mtx_seg.shape[0] == num_segs
    return (
        snps,
        tot_mat,
        ref_mat,
        alt_mat,
        a_mtx,
        b_mtx,
        segs_df,
        a_mtx_seg,
        b_mtx_seg,
        tot_mtx_seg,
    )
