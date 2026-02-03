import os
import numpy as np
import pandas as pd

from scipy.stats import norm


def interp_cM_blocks(
    blocks: pd.DataFrame,
    snp_info: pd.DataFrame,
    genetic_map: pd.DataFrame,
):
    blocks = blocks.copy(deep=True)
    blocks["dist_cM"] = 0.0

    # blocks SNP-start and SNP-end
    hb_pos = snp_info.groupby("bbc_id", sort=False)["POS"].agg(
        snp_start="min", snp_end="max"
    )
    blocks = blocks.join(hb_pos, on="bbc_id")

    genetic_map_chrs = genetic_map.groupby(by="#CHR", sort=False, observed=True)
    for ch, block_idx in blocks.groupby("#CHR", sort=False).groups.items():
        ch_map = genetic_map_chrs.get_group(ch)
        start_cMs = np.interp(
            blocks.loc[block_idx, "snp_start"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )
        end_cMs = np.interp(
            blocks.loc[block_idx, "snp_end"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )

        dist_cM = np.zeros(len(block_idx), dtype=np.float32)
        dist_cM[1:] = start_cMs[1:] - end_cMs[:-1]
        blocks.loc[block_idx, "dist_cM"] = np.maximum(dist_cM, 0.0)
    return blocks["dist_cM"].to_numpy()


def interp_cM_snps(snp_info: pd.DataFrame, genetic_map: pd.DataFrame):
    snp_info["dist_cM"] = 0.0
    genetic_map_chrs = genetic_map.groupby(by="#CHR", sort=False, observed=True)
    for ch, ch_snps in snp_info.groupby(by="#CHR", sort=False):
        ch_maps = genetic_map_chrs.get_group(ch)
        pos_cms = np.interp(
            ch_snps["POS"].to_numpy(),
            ch_maps["POS"].to_numpy(),
            ch_maps["cM"].to_numpy(),
        )

        dist_cM = np.zeros(len(pos_cms), dtype=np.float32)
        dist_cM[1:] = pos_cms[1:] - pos_cms[:-1]
        snp_info.loc[ch_snps.index, "dist_cM"] = np.maximum(dist_cM, 0.0)
    return snp_info["dist_cM"].to_numpy()


def estimate_switchprobs_cM(dist_cms: np.ndarray, nu=1, min_switchprob=1e-6):
    switchprobs = (1 - np.exp(-2 * nu * dist_cms)) / 2.0
    return np.clip(switchprobs, a_min=min_switchprob, a_max=None)


def estimate_switchprobs_PS(blocks: pd.DataFrame, switchprob_ps=0.05):
    switch_bias = 1e-4
    same_block = blocks["PS"] == blocks["PS"].shift(1).fillna(False)
    switchprobs = np.where(
        same_block,
        switchprob_ps,  # within same PS phase block
        0.5 - switch_bias,  # across different PS phase block
    )
    return switchprobs


def binom_2prop_test(
    Y: np.ndarray,
    D: np.ndarray,
    p: np.ndarray,
    alpha: float = 5e-3,
    center: float = 0.5,
    margin: float = 0.05,  # deviation from 0.5
):
    """
    detect phase switch for pairwise SNPs with binomial and sign change from 0.5.
    """
    N = Y.size
    p1 = p[:-1]
    p2 = p[1:]
    d1 = D[:-1]
    d2 = D[1:]
    denom = d1 + d2

    # pooled proportion
    p_pool = np.divide(
        Y[:-1] + Y[1:], denom, out=np.full(N - 1, np.nan), where=denom > 0
    )

    # standard error
    var = p_pool * (1.0 - p_pool) * (1.0 / d1 + 1.0 / d2)
    se = np.sqrt(var)

    # z statistic
    z = np.divide(p2 - p1, se, out=np.full_like(se, np.nan), where=se > 0)

    # critical value for two-sided test
    zcrit = norm.isf(alpha / 2.0)

    # significance per pair
    sig_diff = np.abs(z) > zcrit

    # crossing 0.5 with margin:
    # p1 well below 0.5 & p2 well above 0.5, or vice versa
    low = center - margin
    high = center + margin
    cross_mask = ((p1 < low) & (p2 > high)) | ((p1 > high) & (p2 < low))

    # final decision: must be significant AND cross 0.5 with margin
    pair_diff = sig_diff & cross_mask

    # build S of length N (S[i] is split between i-1 and i)
    S = np.zeros(N, dtype=bool)
    S[1:] = pair_diff & np.isfinite(z)
    return S
