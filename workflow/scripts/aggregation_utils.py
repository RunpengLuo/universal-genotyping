import os
import logging

import numpy as np
import pandas as pd
import numba

from scipy.io import mmread
from scipy.sparse import csr_matrix, hstack, issparse
from scipy.stats import beta

import pyranges as pr
import anndata
import scanpy as sc

from io_utils import *
from postprocess_utils import *


@numba.njit
def _bin_snps_numba(read_counts, min_snp_reads, min_snp_per_block):
    """Greedy adaptive binning over a contiguous read-count array.

    Parameters
    ----------
    read_counts : (N, M) contiguous float64
        Per-SNP read counts for tumor samples.
    min_snp_reads : int
        Minimum total reads per sample for a bin to be complete.
    min_snp_per_block : int
        Minimum number of SNPs per bin.

    Returns
    -------
    bin_ids : (N,) int64
        Relative bin ID for each SNP within this group.
    n_bins : int
        Number of bins created (before last-block adjustment).
    """
    N, M = read_counts.shape
    bin_ids = np.zeros(N, dtype=np.int64)
    if N == 0:
        return bin_ids, 0

    bin_id = 0
    prev_start = 0
    acc = read_counts[0].copy()
    acc_n = 1

    for i in range(1, N):
        # check if accumulator meets both thresholds
        meets_reads = True
        if min_snp_reads > 0:
            for j in range(M):
                if acc[j] < min_snp_reads:
                    meets_reads = False
                    break
        if meets_reads and acc_n >= min_snp_per_block:
            bin_ids[prev_start:i] = bin_id
            bin_id += 1
            prev_start = i
            acc = read_counts[i].copy()
            acc_n = 1
        else:
            for j in range(M):
                acc[j] += read_counts[i, j]
            acc_n += 1

    # last block: keep as own bin if it meets both criteria and isn't the only block
    last_meets_reads = True
    if min_snp_reads > 0:
        for j in range(M):
            if acc[j] < min_snp_reads:
                last_meets_reads = False
                break
    if last_meets_reads and acc_n >= min_snp_per_block and prev_start > 0:
        bin_ids[prev_start:] = bin_id
        bin_id += 1
    else:
        merge_id = bin_id - 1 if bin_id > 0 else 0
        bin_ids[prev_start:] = merge_id

    return bin_ids, bin_id


def _pyranges_assign_chrom(qry, qry_mask, ref_chrom, ref_id, pos_col):
    """PyRanges fallback for assigning positions to overlapping intervals on one chromosome."""
    qry_sub = qry.loc[qry_mask]
    qry_pr = pr.PyRanges(
        chromosomes=qry_sub["#CHR"], starts=qry_sub[pos_col], ends=qry_sub[pos_col] + 1
    )
    qry_pr = qry_pr.insert(pd.Series(data=qry_sub.index.to_numpy(), name="qry_index"))
    ref_pr = pr.PyRanges(
        chromosomes=ref_chrom["#CHR"],
        starts=ref_chrom["START"],
        ends=ref_chrom["END"],
    )
    ref_pr = ref_pr.insert(ref_chrom[ref_id])
    joined = qry_pr.join(ref_pr)
    if len(joined) > 0:
        qry.loc[joined.df["qry_index"], ref_id] = joined.df[ref_id].to_numpy()


def assign_pos_to_range(
    qry: pd.DataFrame,
    ref: pd.DataFrame,
    ref_id="region_id",
    pos_col="POS0",
    nodup=True,
):
    """Assign each query position to the reference interval it falls within.

    Uses ``np.searchsorted`` for non-overlapping intervals (fast path) and
    falls back to PyRanges for chromosomes with overlapping intervals.

    Parameters
    ----------
    qry : pd.DataFrame
        Query DataFrame with ``#CHR`` and *pos_col* columns.
    ref : pd.DataFrame
        Reference intervals with ``#CHR``, ``START``, ``END``, and *ref_id*.
    ref_id : str
        Column name for the reference interval identifier.
    pos_col : str
        Column in *qry* containing 0-based positions.
    nodup : bool
        If True, each query is assigned to at most one reference interval;
        if False, returns all overlapping (query, ref) pairs.

    Returns
    -------
    pd.DataFrame
        *qry* with *ref_id* column added (if ``nodup=True``), or a hits
        DataFrame (if ``nodup=False``).
    """
    if not nodup:
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
        hits = (
            qry_pr.join(ref_pr)
            .df[["Chromosome", "Start", "qry_index", ref_id]]
            .rename(columns={"Chromosome": "#CHR", "Start": "POS0"})
        )
        hits["POS"] = hits["POS0"] + 1
        return hits

    # nodup=True: assign each query position to at most one interval
    qry[ref_id] = pd.NA
    for chrom in ref["#CHR"].unique():
        qry_mask = (qry["#CHR"] == chrom).to_numpy()
        if not qry_mask.any():
            continue

        ref_chrom = ref.loc[ref["#CHR"] == chrom].sort_values("START")
        starts = ref_chrom["START"].to_numpy()
        ends = ref_chrom["END"].to_numpy()
        ids = ref_chrom[ref_id].to_numpy()
        positions = qry.loc[qry_mask, pos_col].to_numpy()

        has_overlap = len(starts) > 1 and np.any(starts[1:] < ends[:-1])

        if not has_overlap:
            # Fast path: searchsorted for non-overlapping intervals
            # 0-based half-open: START <= pos < END
            idx = np.searchsorted(starts, positions, side="right") - 1
            safe_idx = idx.clip(min=0)
            valid = (idx >= 0) & (positions < ends[safe_idx])
            qry_indices = qry.index[qry_mask]
            qry.loc[qry_indices[valid], ref_id] = ids[idx[valid]]
        else:
            # Fallback: PyRanges join for overlapping intervals on this chromosome
            _pyranges_assign_chrom(qry, qry_mask, ref_chrom, ref_id, pos_col)

    return qry


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
    logging.info(f"adaptive binning, colname={colname}, min_snp_reads={min_snp_reads}, min_snp_per_block={min_snp_per_block}")
    logging.info(f"#SNP groups={len(snp_grps)}, grouper: {grp_cols}")
    for _, grp_snps in snp_grps:
        grp_idxs = grp_snps.index.to_numpy()

        # Extract the group's read-count submatrix as a contiguous dense copy
        grp_mtx = tot_mtx[grp_idxs]
        if issparse(grp_mtx):
            grp_reads = np.ascontiguousarray(grp_mtx[:, tumor_sidx:].toarray(), dtype=np.float64)
        else:
            grp_reads = np.ascontiguousarray(grp_mtx[:, tumor_sidx:], dtype=np.float64).copy()

        # Call the Numba kernel for greedy binning
        local_bin_ids, n_bins = _bin_snps_numba(grp_reads, min_snp_reads, min_snp_per_block)

        # Offset local bin IDs to global bin IDs
        local_bin_ids += bin_id
        snps.loc[grp_idxs, colname] = local_bin_ids
        bin_id += max(n_bins, 1)

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
    for feat_col in ("feature_id", "feature_type"):
        if feat_col in snps.columns:
            snp_bins[feat_col] = snp_grps[feat_col].agg(
                lambda x: ";".join(dict.fromkeys(str(v) for v in x if pd.notna(v)))
            )
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

    _rename = {k: v for k, v in {"#CHR": "Chromosome", "START": "Start", "END": "End"}.items()
               if k in qry.columns and v not in qry.columns}
    qry_gr = pr.PyRanges(qry.rename(columns=_rename))
    _rename = {k: v for k, v in {"#CHR": "Chromosome", "START": "Start", "END": "End"}.items()
               if k in ref.columns and v not in ref.columns}
    ref_gr = pr.PyRanges(ref.rename(columns=_rename))

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
