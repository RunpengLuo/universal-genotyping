import os
import logging

import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy.sparse import csr_matrix, hstack, issparse
from scipy.stats import beta

import pyranges as pr

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import *


##################################################
def canon_mat_one_replicate(
    parent_keys: pd.Index,
    vcf_file: str,
    tot_mtx_file: str,
    ad_mtx_file: str,
    ncells: int,
):
    """
    Map SNP by barcode DP/AD sparse mats to same SNP position index defined by parent SNP file.
    """
    M = len(parent_keys)
    child_snps = read_VCF(vcf_file, addkey=True)
    m = len(child_snps)

    tot_mtx: csr_matrix = mmread(tot_mtx_file).tocsr()
    ad_mtx: csr_matrix = mmread(ad_mtx_file).tocsr()
    assert tot_mtx.shape == ad_mtx.shape
    assert tot_mtx.shape == (m, ncells)

    # reorder mats to preserve same order as child_snps (sorted when read_VCF)
    raw_snp_ids = child_snps["RAW_SNP_IDX"].to_numpy()
    tot_mtx = tot_mtx[raw_snp_ids, :]
    ad_mtx = ad_mtx[raw_snp_ids, :]

    dup_mask = child_snps["KEY"].duplicated(keep=False).to_numpy()
    n_dup_rows = int(dup_mask.sum())
    if n_dup_rows > 0:
        n_dup_keys = int(child_snps.loc[dup_mask, "KEY"].nunique())
        logging.warning(
            f"#found duplicated SNP #CHR/POS, #rows={n_dup_rows}/{m} #SNPs={n_dup_keys}, drop all"
        )
        child_snps = child_snps.loc[~dup_mask].reset_index(drop=True)
        tot_mtx = tot_mtx[~dup_mask, :]
        ad_mtx = ad_mtx[~dup_mask, :]

    m = len(child_snps)
    child_keys = pd.Index(child_snps["KEY"])

    # for each parent SNP, find its index in child SNP, -1 if not presented.
    # e.g., [0, -1, 1, 2, -1], M=5, m=3
    child_loc = child_keys.get_indexer(parent_keys)
    # e.g., [1, 0, 1, 1, 0]
    present = child_loc >= 0
    logging.info(f"located SNPs in parent={present.sum()}/{M}")

    # parent-row indices that exist in child slice
    # e.g., [0, 2, 3], located SNP index in parent [0, M).
    common_pidx = np.flatnonzero(present)
    # e.g., [0, 1, 2], located SNP index in child [0, m)
    common_cidx = child_loc[present]

    # expanded DP mtx
    tot_present = tot_mtx[common_cidx, :].tocoo()
    tot_canon = csr_matrix(
        (tot_present.data, (common_pidx[tot_present.row], tot_present.col)),
        shape=(M, ncells),
        dtype=tot_mtx.dtype,
    )

    # expanded AD mtx
    ad_present = ad_mtx[common_cidx, :].tocoo()
    ad_canon = csr_matrix(
        (ad_present.data, (common_pidx[ad_present.row], ad_present.col)),
        shape=(M, ncells),
        dtype=ad_mtx.dtype,
    )

    return tot_canon, ad_canon


def merge_mats(tot_list: list, ad_list: list):
    """Horizontally stack per-replicate total and alt-count matrices.

    Derives the REF matrix as ``TOT - ALT``.

    Parameters
    ----------
    tot_list : list of csr_matrix
        Per-replicate total depth matrices (SNPs x cells).
    ad_list : list of csr_matrix
        Per-replicate alt-allele count matrices (SNPs x cells).

    Returns
    -------
    tuple[csr_matrix, csr_matrix, csr_matrix]
        ``(tot_mtx, ref_mtx, alt_mtx)`` concatenated across replicates.
    """
    if len(tot_list) == 1:
        return tot_list[0], tot_list[0] - ad_list[0], ad_list[0]
    tot_mtx = hstack(tot_list, format="csr")
    alt_mtx = hstack(ad_list, format="csr")
    ref_mtx = tot_mtx - alt_mtx
    return tot_mtx, ref_mtx, alt_mtx


def hstack_mtx_list(mtx_list: list):
    """Horizontally stack a list of sparse matrices into a single CSR matrix.

    Parameters
    ----------
    mtx_list : list of csr_matrix
        Matrices with the same number of rows.

    Returns
    -------
    csr_matrix
        Concatenated matrix.
    """
    if len(mtx_list) == 1:
        return mtx_list[0]
    return hstack(mtx_list, format="csr")


def apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, phases):
    """Apply per-SNP phase labels to produce phased A/B allele count matrices.

    Parameters
    ----------
    tot_mtx : sparse or ndarray
        Total depth matrix (SNPs x cells/samples).
    ref_mtx : sparse or ndarray
        Reference allele count matrix.
    alt_mtx : sparse or ndarray
        Alternate allele count matrix.
    phases : np.ndarray
        Per-SNP phase labels (0 or 1).

    Returns
    -------
    tuple
        ``(a_mtx, b_mtx)`` — phased allele count matrices.
    """
    p = phases[:, None]
    if issparse(ref_mtx):
        b_mtx = ref_mtx.multiply(p) + alt_mtx.multiply(1 - p)
        b_mtx.data = np.rint(b_mtx.data).astype(np.int32)
    else:
        b_mtx = ref_mtx * p + alt_mtx * (1 - p)
        b_mtx = np.round(b_mtx).astype(np.int32)
    a_mtx = tot_mtx - b_mtx
    return a_mtx, b_mtx


def compute_af_per_sample(tot_mtx, b_mtx, i: int):
    """Compute per-SNP allele frequency for a single sample column.

    Parameters
    ----------
    tot_mtx : sparse or ndarray
        Total depth matrix (SNPs x samples).
    b_mtx : sparse or ndarray
        B-allele count matrix.
    i : int
        Sample column index.

    Returns
    -------
    np.ndarray
        Allele frequency per SNP; ``NaN`` where depth is zero.
    """
    tot_col = tot_mtx[:, i]
    b_col = b_mtx[:, i]

    den = (
        tot_col.toarray().ravel() if issparse(tot_col) else np.asarray(tot_col).ravel()
    )
    num = b_col.toarray().ravel() if issparse(b_col) else np.asarray(b_col).ravel()

    out = np.full_like(den, np.nan, dtype=np.float32)
    return np.divide(num, den, out=out, where=(den > 0))


def compute_af_pseudobulk(tot_mtx, b_mtx):
    """Compute per-SNP allele frequency across all cells (pseudobulk sum).

    Parameters
    ----------
    tot_mtx : sparse or ndarray
        Total depth matrix (SNPs x cells).
    b_mtx : sparse or ndarray
        B-allele count matrix.

    Returns
    -------
    np.ndarray
        Pseudobulk allele frequency per SNP; ``NaN`` where total depth is zero.
    """
    if issparse(tot_mtx):
        den = np.asarray(tot_mtx.sum(axis=1)).ravel()
    else:
        den = tot_mtx.sum(axis=1)

    if issparse(b_mtx):
        num = np.asarray(b_mtx.sum(axis=1)).ravel()
    else:
        num = b_mtx.sum(axis=1)

    out = np.full(den.shape[0], np.nan, dtype=np.float32)
    return np.divide(num, den, out=out, where=(den > 0))


##################################################
def get_mask_by_region(snps: pd.DataFrame, regions: pd.DataFrame) -> np.ndarray:
    """
    Return a boolean mask (len == len(snps)) indicating whether each SNP (CHR, POS)
    overlaps any interval in a BED-like file (Chromosome, Start, End), using 0-based
    half-open intervals [Start, End).

    Assumes SNP POS is 1-based. Converts each SNP to an interval [POS-1, POS).
    """
    snp_positions = snps[["#CHR", "POS"]].copy()
    snp_positions["Start"] = snp_positions["POS"].astype(np.int64) - 1
    snp_positions["End"] = snp_positions["POS"].astype(np.int64)

    pr_snps = pr.PyRanges(snp_positions.rename(columns={"#CHR": "Chromosome"}))
    pr_regions = pr.PyRanges(regions)

    overlapping = pr_snps.overlap(pr_regions).df.rename(columns={"Chromosome": "#CHR"})

    mask = (
        snp_positions[["#CHR", "POS"]]
        .merge(
            overlapping[["#CHR", "POS"]], on=["#CHR", "POS"], how="left", indicator=True
        )["_merge"]
        .eq("both")
        .to_numpy(dtype=bool)
    )
    return mask


def get_mask_by_depth(snps: pd.DataFrame, tot_mtx: csr_matrix, min_dp=1):
    """Return a boolean mask keeping SNPs where every sample meets the minimum depth.

    Parameters
    ----------
    snps : pd.DataFrame
        SNP info DataFrame (used only for logging).
    tot_mtx : csr_matrix
        Total depth matrix (SNPs x samples).
    min_dp : int
        Minimum depth threshold per sample.

    Returns
    -------
    np.ndarray
        Boolean mask of length ``len(snps)``.
    """
    mask = np.all(tot_mtx >= min_dp, axis=1)
    logging.info(
        f"filter by depth, min_dp={min_dp}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask


def get_mask_by_depth_pseudobulk(snps: pd.DataFrame, tot_mtx: csr_matrix, min_dp=1):
    """Return a boolean mask keeping SNPs whose pseudobulk (row-summed) depth meets the threshold.

    Parameters
    ----------
    snps : pd.DataFrame
        SNP info DataFrame (used only for logging).
    tot_mtx : csr_matrix
        Total depth matrix (SNPs x cells).
    min_dp : int
        Minimum pseudobulk depth threshold.

    Returns
    -------
    np.ndarray
        Boolean mask of length ``len(snps)``.
    """
    mask = np.asarray(tot_mtx.sum(axis=1)).ravel() >= min_dp
    logging.info(
        f"filter by depth, min_dp={min_dp}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask


def get_mask_by_het_balanced(
    snps: pd.DataFrame,
    ref_mtx: csr_matrix,
    alt_mtx: csr_matrix,
    gamma: float,
    normal_idx=0,
):
    """
    mask SNPs if normal sample failed beta-posterior credible interval test with beta(1, 1) prior.
    """
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    q = np.array([p_lower, p_upper])
    het_cred_ints = beta.ppf(
        q[None, :],
        ref_mtx[:, normal_idx][:, None] + 1,
        alt_mtx[:, normal_idx][:, None] + 1,
    )
    mask = (het_cred_ints[:, 0] <= 0.5) & (0.5 <= het_cred_ints[:, 1])
    logging.info(
        f"filter by balanced Het-SNPs on normal sample, gamma={gamma}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask


##################################################
def subset_baf(
    baf_df: pd.DataFrame, ch: str, start: int, end: int, is_last_block=False
):
    """Slice a BAF DataFrame to a chromosomal interval ``[start, end)``.

    For the last block, the interval is closed on the right: ``[start, end]``.

    Parameters
    ----------
    baf_df : pd.DataFrame
        DataFrame with ``#CHR`` and ``POS`` columns (or ``POS`` as index).
    ch : str or None
        Chromosome to filter on; if None, no chromosome filter is applied.
    start, end : int
        Genomic position boundaries.
    is_last_block : bool
        If True, use a closed right boundary.

    Returns
    -------
    pd.DataFrame
        Filtered subset.
    """
    if ch != None:
        baf_ch = baf_df[baf_df["#CHR"] == ch]
    else:
        baf_ch = baf_df
    if baf_ch.index.name == "POS":
        pos = baf_ch.index
    else:
        pos = baf_ch["POS"]
    if is_last_block:
        return baf_ch[(pos >= start) & (pos <= end)]
    else:
        return baf_ch[(pos >= start) & (pos < end)]


def assign_snp_bounderies(
    snps: pd.DataFrame, regions: pd.DataFrame, colname="region_id"
):
    """
    divide regions into [START, END) subregions, each subregion has one SNP.
    If a SNP is out-of-region, its START and END will be 0 and region_id will be 0.
    """
    snps["START"] = 0
    snps["END"] = 0

    snps[colname] = 0
    region_id = 0

    chroms = snps["#CHR"].unique().tolist()
    region_grps_ch = regions.groupby(by="#CHR", sort=False)
    for chrom in chroms:
        regions_ch = region_grps_ch.get_group(chrom)
        for _, region in regions_ch.iterrows():
            reg_start, reg_end = region["START"], region["END"]
            reg_snps = subset_baf(snps, chrom, reg_start, reg_end)
            if len(reg_snps) == 0:
                continue
            reg_snp_positions = reg_snps["POS0"].to_numpy()
            reg_snp_indices = reg_snps.index.to_numpy()

            # annotate region ID
            snps.loc[reg_snp_indices, colname] = region_id
            region_id += 1

            # build SNP bounderies
            if len(reg_snps) == 1:
                snps.loc[reg_snp_indices, "START"] = reg_start
                snps.loc[reg_snp_indices, "END"] = reg_end
            else:
                reg_bounderies = np.ceil(
                    np.vstack([reg_snp_positions[:-1], reg_snp_positions[1:]]).mean(
                        axis=0
                    )
                ).astype(np.uint32)
                reg_bounderies = np.concatenate(
                    [[reg_start], reg_bounderies, [reg_end]]
                )
                snps.loc[reg_snp_indices, "START"] = reg_bounderies[:-1]
                snps.loc[reg_snp_indices, "END"] = reg_bounderies[1:]

    snps["BLOCKSIZE"] = snps["END"] - snps["START"]
    return snps


##################################################
def plot_allele_freqs(
    pos_df,
    rep_ids,
    tot_mtx,
    b_mtx,
    genome_size,
    plot_dir,
    apply_pseudobulk=False,
    allele="ref",
    unit="SNP",
    suffix="",
):
    """Generate genome-wide allele-frequency scatter plots.

    Produces either one plot per sample or a single pseudobulk plot, saved as
    multi-page PDFs (one page per chromosome).

    Parameters
    ----------
    pos_df : pd.DataFrame
        SNP/bin position DataFrame with ``#CHR`` and ``POS`` (or ``START``/``END``).
    rep_ids : list[str]
        Replicate identifiers.
    tot_mtx, b_mtx : sparse or ndarray
        Total depth and B-allele count matrices.
    genome_size : str
        Path to chromosome sizes file.
    plot_dir : str
        Output directory for PDF plots.
    apply_pseudobulk : bool
        If True, plot a single pseudobulk AF; otherwise plot per-sample.
    allele : str
        Allele label for filenames (e.g., ``"ref"``, ``"B"``).
    unit : str
        Feature unit label (e.g., ``"SNP"``, ``"bin"``).
    suffix : str
        Optional filename suffix.
    """
    logging.info(
        f"QC analysis - plot {allele}-{unit} allele frequency, {unit}={allele}, apply_pseudobulk={apply_pseudobulk}"
    )
    if apply_pseudobulk:
        af = compute_af_pseudobulk(tot_mtx, b_mtx)
        plot_file = os.path.join(plot_dir, f"af_{allele}_{unit}.pseudobulk{suffix}.pdf")
        plot_1d_sample(pos_df, af, genome_size, plot_file, unit=unit, val_type="AF")
    else:
        _tot_mtx = tot_mtx.tocsc() if issparse(tot_mtx) else tot_mtx
        _b_mtx = b_mtx.tocsc() if issparse(b_mtx) else b_mtx
        for i, rep_id in enumerate(rep_ids):
            plot_file = os.path.join(
                plot_dir, f"af_{allele}_{unit}.{rep_id}{suffix}.pdf"
            )
            af = compute_af_per_sample(_tot_mtx, _b_mtx, i)
            plot_1d_sample(pos_df, af, genome_size, plot_file, unit=unit, val_type="AF")
    return


def plot_1d_sample(
    pos_df: pd.DataFrame,
    val: np.ndarray,
    genome_size: str,
    out_file: str,
    unit="SNP",
    val_type="BAF",
    s=4,
    dpi=150,
    alpha=0.6,
    figsize=(40, 3),
    min_ylim=0.0,
    max_ylim=1.0,
):
    """
    plot any features like AF, RDR, etc., in 1d chromosome scatter plot.
    """
    logging.info(f"chromosome wide {unit}-level {val_type} plot, out_file={out_file}")
    chrom_sizes = get_chr_sizes(genome_size)

    ch = pos_df["#CHR"].to_numpy()
    if "POS" in pos_df.columns:
        pos = pos_df["POS"].to_numpy()
    else:
        # use midpoints for [START, END) intervals
        pos = ((pos_df["START"].to_numpy() + pos_df["END"].to_numpy()) // 2).astype(
            np.int64
        )
    num_bins = len(pos_df)
    # find contiguous chromosome blocks
    change = np.flatnonzero(ch[1:] != ch[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, num_bins]
    chroms = ch[starts]

    pdf_fd = PdfPages(out_file)
    for chrom, lo, hi in zip(chroms, starts, ends):
        chr_end = chrom_sizes.get(chrom)
        if chr_end is None:
            logging.warning(f"{chrom}: not found in {genome_size}")
            continue

        x = pos[lo:hi]
        y = val[lo:hi]
        m = np.isfinite(y)

        n_plot = int(m.sum())
        if n_plot == 0:
            logging.warning(f"{chrom}: all {val_type} values are non-finite")
            continue
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.scatter(x[m], y[m], s=s, alpha=alpha, rasterized=True)
        if val_type in ["AF", "BAF"]:
            ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
            ax.set_ylim(-0.05, 1.05)
        else:
            ax.set_ylim(min_ylim, max_ylim)

        ax.set_xlabel(val_type)
        ax.set_xlim(0, chr_end)
        ax.grid(alpha=0.2)
        ax.set_title(f"{val_type} plot - {chrom}")
        pdf_fd.savefig(fig, dpi=dpi)
        plt.close(fig)
    pdf_fd.close()
    return
