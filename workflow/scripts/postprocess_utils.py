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

from utils import read_VCF, get_chr_sizes


##################################################
def canon_mat_one_replicate(
    parent_keys: pd.Index,
    vcf_file: str,
    dp_mat_file: str,
    ad_mat_file: str,
    ncells: int,
):
    """
    Map SNP by barcode DP/AD sparse mats to same SNP position index defined by parent SNP file.
    """
    M = len(parent_keys)
    child_snps = read_VCF(vcf_file, addkey=True)
    m = len(child_snps)
    child_keys = pd.Index(child_snps["KEY"])
    assert not child_keys.duplicated().any()

    # for each parent SNP, find its index in child SNP, -1 if not presented.
    # e.g., [0, -1, 1, 2, -1], M=5, m=3
    child_loc = child_keys.get_indexer(parent_keys)
    # e.g., [1, 0, 1, 1, 0]
    present = child_loc >= 0
    logging.info(f"located SNPs in parent={present.sum()}/{M}")

    dp_mtx: csr_matrix = mmread(dp_mat_file).tocsr()
    ad_mtx: csr_matrix = mmread(ad_mat_file).tocsr()

    assert dp_mtx.shape == ad_mtx.shape
    assert dp_mtx.shape == (m, ncells)

    # parent-row indices that exist in child slice
    # e.g., [0, 2, 3], located SNP index in parent [0, M).
    common_pidx = np.flatnonzero(present)
    # e.g., [0, 1, 2], located SNP index in child [0, m)
    common_cidx = child_loc[present]

    # expanded DP mtx
    dp_present = dp_mtx[common_cidx, :].tocoo()
    dp_canon = csr_matrix(
        (dp_present.data, (common_pidx[dp_present.row], dp_present.col)),
        shape=(M, ncells),
        dtype=dp_mtx.dtype,
    )

    # expanded AD mtx
    ad_present = ad_mtx[common_cidx, :].tocoo()
    ad_canon = csr_matrix(
        (ad_present.data, (common_pidx[ad_present.row], ad_present.col)),
        shape=(M, ncells),
        dtype=ad_mtx.dtype,
    )

    return dp_canon, ad_canon


def merge_mats(dp_list: list, ad_list: list):
    if len(dp_list) == 1:
        return dp_list[0], dp_list[0] - ad_list[0], ad_list[0]
    dp_mtx = hstack(dp_list, format="csr")
    alt_mtx = hstack(ad_list, format="csr")
    ref_mtx = dp_mtx - alt_mtx
    return dp_mtx, ref_mtx, alt_mtx


def apply_phase_to_mat(dp_mtx, ref_mtx, alt_mtx, phases):
    p = phases[:, None]
    if issparse(ref_mtx):
        b_mtx = ref_mtx.multiply(p) + alt_mtx.multiply(1 - p)
        a_mtx = dp_mtx - b_mtx
    else:
        b_mtx = ref_mtx * p + alt_mtx * (1 - p)
        a_mtx = dp_mtx - b_mtx
    return a_mtx, b_mtx


def compute_af(dp_mtx, b_mtx, apply_pseudobulk=False):
    dp = dp_mtx.toarray() if issparse(dp_mtx) else dp_mtx
    b = b_mtx.toarray() if issparse(b_mtx) else b_mtx
    nsnps = dp.shape[0]
    if apply_pseudobulk:
        af = np.divide(
            np.sum(b, axis=1),
            np.sum(dp, axis=1),
            out=np.full(nsnps, np.nan, dtype=np.float32),
        )
    else:
        af = np.divide(b, dp, out=np.full(dp.shape, np.nan, dtype=np.float32))
    return af


##################################################
def get_mask_by_region(snps: pd.DataFrame, region_bed_file: str) -> np.ndarray:
    """
    Return a boolean mask (len == len(snps)) indicating whether each SNP (CHR, POS)
    overlaps any interval in a BED-like file (Chromosome, Start, End), using 0-based
    half-open intervals [Start, End).

    Assumes SNP POS is 1-based. Converts each SNP to an interval [POS-1, POS).
    """
    regions = pd.read_table(
        region_bed_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["Chromosome", "Start", "End"],
        dtype={0: "string"},
    )
    if not str(regions["Chromosome"].iloc[0]).startswith("chr"):
        if str(snps["#CHR"].iloc[0]).startswith("chr"):
            regions["Chromosome"] = "chr" + regions["Chromosome"].astype(str)

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

    logging.info(f"filter by regions, #passed SNPs={int(mask.sum())}/{len(snps)}")
    return mask


def get_mask_by_depth(snps: pd.DataFrame, dp_mtx: csr_matrix, min_dp=1):
    mask = np.asarray(dp_mtx.sum(axis=1)).ravel() >= min_dp
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
    assume first sample is normal.
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
def plot_snps_allele_freqs(
    snps,
    rep_ids,
    dp_mtx,
    b_mtx,
    genome_file,
    plot_dir,
    apply_pseudobulk=False,
    allele="ref",
):
    af = compute_af(dp_mtx, b_mtx, apply_pseudobulk)
    if apply_pseudobulk:
        plot_file = os.path.join(plot_dir, f"af_{allele}_allele.pseudobulk.pdf")
        plot_snps_allele_freqs_sample(snps, af, genome_file, plot_file)
    else:
        for i, rep_id in enumerate(rep_ids):
            plot_file = os.path.join(
                plot_dir, f"af_{allele}_allele.pseudobulk.{rep_id}.pdf"
            )
            plot_snps_allele_freqs_sample(snps, af[:, i], genome_file, plot_file)
    return


def plot_snps_allele_freqs_sample(
    snps: pd.DataFrame,
    af: np.ndarray,
    genome_file: str,
    out_file: str,
    s=4,
    dpi=150,
):
    logging.info(f"plot 1D per-SNP B-allele frequency, out_file={out_file}")
    chrom_sizes = get_chr_sizes(genome_file)

    ch = snps["#CHR"].to_numpy()
    pos = snps["POS"].to_numpy()
    # find contiguous chromosome blocks
    change = np.flatnonzero(ch[1:] != ch[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(snps)]
    chroms = ch[starts]

    pdf_fd = PdfPages(out_file)
    for chrom, lo, hi in zip(chroms, starts, ends):
        chr_end = chrom_sizes.get(chrom)
        if chr_end is None:
            logging.warning(f"{chrom}: not found in {genome_file}")
            continue

        x = pos[lo:hi]
        y = af[lo:hi]
        m = np.isfinite(y)

        n_all = hi - lo
        n_plot = int(m.sum())
        if n_plot == 0:
            logging.warning(f"{chrom}: all SNPs have non-finite AF")
            continue

        logging.info(f"plot {chrom} with #SNP={n_plot}/{n_all}")
        fig, ax = plt.subplots(1, 1, figsize=(40, 3))
        ax.scatter(x[m], y[m], s=s, alpha=0.6, rasterized=True)
        ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlim(0, chr_end)
        ax.grid(alpha=0.2)
        ax.set_title(f"allele-frequency plot - {chrom}")
        pdf_fd.savefig(fig, dpi=dpi)
        plt.close(fig)
    pdf_fd.close()
    return
