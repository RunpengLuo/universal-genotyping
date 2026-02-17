import os
import logging

import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy.sparse import csr_matrix, hstack, issparse
from scipy.stats import beta

import pyranges as pr
import anndata
import scanpy as sc

from io_utils import *
from postprocess_utils import *


def assign_pos_to_range(
    qry: pd.DataFrame,
    ref: pd.DataFrame,
    ref_id="region_id",
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
    snp_df: pd.DataFrame, region_df: pd.DataFrame, assay_type: str, region_id="BIN_ID"
):
    """
    region_df must be 0-indexed non-overlapping intervals [s, t) in standard BED format.
    """
    logging.info(f"#{assay_type}-SNP (raw)={len(snp_df)}")
    snp_df = assign_pos_to_range(snp_df, region_df, ref_id=region_id, pos_col="POS0")
    isna_snp_df = snp_df[region_id].isna()
    logging.info(
        f"#{assay_type}: #SNPS outside any region={np.sum(isna_snp_df) / len(snp_df):.3%}"
    )
    snp_df.dropna(subset=region_id, inplace=True)
    snp_df[region_id] = snp_df[region_id].astype(region_df[region_id].dtype)
    logging.info(f"#{assay_type}-SNP (remain)={len(snp_df)}")

    counts = snp_df[region_id].value_counts()
    region_df[f"#SNPS"] = region_df[region_id].map(counts).fillna(0).astype(int)
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
    Adaptive binning over pre-defined chromosome regions.
    Each bin satisfies both MSR (min total reads) and MSPB (min SNPs per block).
    The last block in each group is kept as its own bin if it meets both criteria,
    otherwise it is merged into the previous bin.
    In-place adds <colname> column to snps to indicate bin ids.
    """
    bin_id = 0
    snps[colname] = 0
    snp_grps = snps.groupby(by=grp_cols, sort=False)
    logging.info(f"adaptive binning, colname={colname}")
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

        # fill last block: keep as its own bin if it meets both criteria,
        # otherwise merge backward into the previous bin.
        if (
            np.all(acc_read_count >= min_snp_reads)
            and acc_num_snp >= min_snp_per_block
            and prev_start > 0
        ):
            bin_ids[prev_start:] = bin_id
            bin_id += 1
        else:
            bin_ids[prev_start:] = max(bin_id - 1, bin_id0)
        snps.loc[grp_idxs, colname] = bin_ids

        bin_id = max(bin_id, bin_id0 + 1)

    pos_dict = {
        "#CHR": ("#CHR", "first"),
        "START": ("START", "min"),
        "END": ("END", "max"),
        # effective first/last SNP positions.
        "START0": ("POS0", "min"),
        "END0": ("POS", "max"),
    }
    # only include switchprobs if present (not guaranteed when called externally)
    if "switchprobs" in snps.columns:
        pos_dict["switchprobs"] = ("switchprobs", "first")
    for grp_col in grp_cols:
        pos_dict[grp_col] = (grp_col, "first")

    snp_grps = snps.groupby(by=colname, sort=False, as_index=True)
    snp_bins = snp_grps.agg(**pos_dict)
    snp_bins.loc[:, "#SNPS"] = snp_grps.size()  # align by bin_id index, not position
    snp_bins.loc[:, "BLOCKSIZE"] = snp_bins["END"] - snp_bins["START"]
    snp_bins[colname] = snp_bins.index

    bin_sizes = snp_bins["#SNPS"].to_numpy()
    block_sizes = snp_bins["BLOCKSIZE"].to_numpy()
    logging.info("adaptive binning summary")
    logging.info(f"#SNPs={len(snps)}")
    logging.info(f"#bins={len(snp_bins)}")
    logging.info(
        "snps per bin: min=%.0f  median=%.0f  max=%.0f",
        float(bin_sizes.min()),
        float(np.median(bin_sizes)),
        float(bin_sizes.max()),
    )
    logging.info(
        "blocksize per bin: min=%.0f  median=%.0f  max=%.0f",
        float(block_sizes.min()),
        float(np.median(block_sizes)),
        float(block_sizes.max()),
    )

    return snp_bins


def matrix_segmentation(X, bin_ids, K):
    """
    N: #features
    M: #samples
    X: (N, M) sparse or dense   [snp-by-sample]
    bin_ids: (N,) ints in [0..K-1]  (assign each SNP to a bin, bin ids are non-decreasing)
    return:
      - sparse in  -> (K, M) csr_matrix
      - dense in   -> (K, M) ndarray
    """
    X = X.tocsr() if issparse(X) else np.asarray(X)

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
    return X_bin


def assign_largest_overlap(
    qry: pd.DataFrame, ref: pd.DataFrame, qry_id: str, ref_id: str
) -> pd.DataFrame:
    """
    for each row in qry, assign the ID of the overlapping interval in ref
    that has largest overlap length
    """

    qry_gr = pr.PyRanges(
        qry.rename(columns={"#CHR": "Chromosome", "START": "Start", "END": "End"})
    )
    ref_gr = pr.PyRanges(
        ref.rename(columns={"#CHR": "Chromosome", "START": "Start", "END": "End"})
    )

    joined = qry_gr.join(ref_gr, suffix="_REF").as_df()
    joined["overlap_len"] = (
        np.minimum(joined["End"], joined["End_REF"])
        - np.maximum(joined["Start"], joined["Start_REF"])
    ).clip(lower=0)

    best = joined.loc[joined.groupby(joined.index)["overlap_len"].idxmax()].copy()

    best = best[[qry_id, ref_id, "overlap_len"]]
    qry = qry.merge(best, on=[qry_id], how="left", sort=False)
    max_ovlp_idx = qry.groupby(by=qry_id, sort=False)["overlap_len"].apply(
        lambda s: s.idxmax() if s.notna().any() else s.index[0]
    )
    qry_out = qry.loc[max_ovlp_idx].sort_values(by=qry_id).reset_index(drop=True)
    qry_out = qry_out.drop(columns="overlap_len")
    return qry_out


def feature_to_blocks(
    adata: sc.AnnData,
    blocks: pd.DataFrame,
    assay_type: str,
    feature_idx="feature_idx",
    block_idx="region_id",
    drop_cols=True,
):
    """
    filter features not in blocks, likely masked regions include centromeres
    """
    logging.info(f"assign {assay_type} features to blocks, {feature_idx}-{block_idx}")
    adata.var[feature_idx] = np.arange(len(adata.var))

    feature_df = adata.var.reset_index(drop=True)
    logging.info(f"#{assay_type}-features (raw)={len(feature_df)}")

    feature_df = assign_largest_overlap(feature_df, blocks, feature_idx, block_idx)
    isna_features = feature_df[block_idx].isna()
    logging.info(
        f"#{assay_type} feature outside any blocks={np.sum(isna_features) / len(feature_df):.3%}"
    )
    feature_df.dropna(subset=block_idx, inplace=True)
    feature_df[block_idx] = feature_df[block_idx].astype(np.int32)
    logging.info(f"#{assay_type} feature (remain)={len(feature_df)}")

    ##################################################
    # map HB tag to adata.var, filter any out-of-range features
    adata.var = (
        adata.var.reset_index(drop=False)
        .merge(
            right=feature_df[[feature_idx, block_idx]],
            on=feature_idx,
            how="left",
        )
        .set_index("index")
    )
    adata = adata[:, adata.var[block_idx].notna()].copy()
    if drop_cols:
        adata.var.drop(columns=[feature_idx, block_idx], inplace=True)
    else:
        adata.var[block_idx] = adata.var[block_idx].astype(feature_df[block_idx].dtype)
    return adata
