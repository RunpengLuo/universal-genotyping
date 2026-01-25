import logging

import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy.sparse import csr_matrix, hstack, issparse
from scipy.stats import beta

import pyranges as pr

from utils import read_VCF


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

##################################################
def get_mask_by_region(
    snps: pd.DataFrame,
    region_bed_file: str,
):
    regions = pd.read_table(
        region_bed_file,
        sep="\t",
        header=None,
        usecols=range(3),
        names=["Chromosome", "Start", "End"],
    )
    snp_positions = snps[["#CHR", "POS"]].reset_index(drop=True)
    # 0-indexed [Start, End) region
    snp_positions["Start"] = snp_positions["POS"] - 1
    snp_positions["End"] = snp_positions["POS"]

    pr_snps = pr.PyRanges(snp_positions.rename(columns={"#CHR": "Chromosome"}))
    pr_regions = pr.PyRanges(regions)

    overlapping_snps = pr_snps.overlap(pr_regions).df
    overlapping_snps = overlapping_snps.rename(columns={"Chromosome": "#CHR"})
    mask = (
        pd.merge(
            left=snp_positions,
            right=overlapping_snps,
            on=["#CHR", "POS"],
            how="left",
            sort=False,
        )["Start"]
        .notna()
        .to_numpy()
    )
    logging.info(f"filter by regions, #passed SNPs={np.sum(mask)}/{len(snps)}")
    return mask


def get_mask_by_depth(snps: pd.DataFrame, dp_mtx: csr_matrix, min_dp=1):
    mask = np.asarray(dp_mtx.sum(axis=1)).ravel() >= min_dp
    logging.info(
        f"filter by depth, min_dp={min_dp}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask


def get_mask_by_het_balanced(
    snps: pd.DataFrame, ref_mtx: csr_matrix, alt_mtx: csr_matrix, gamma: float, normal_idx=0
):
    """
    mask SNPs if normal sample failed beta-posterior credible interval test with beta(1, 1) prior.
    assume first sample is normal.
    """
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    q = np.array([p_lower, p_upper])
    het_cred_ints = beta.ppf(
        q[None, :], ref_mtx[:, normal_idx][:, None] + 1, alt_mtx[:, normal_idx][:, None] + 1
    )
    mask = (het_cred_ints[:, 0] <= 0.5) & (0.5 <= het_cred_ints[:, 1])
    logging.info(
        f"filter by balanced Het-SNPs on normal sample, gamma={gamma}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask
