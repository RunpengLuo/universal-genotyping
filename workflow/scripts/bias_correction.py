import os
import logging

import numpy as np
import pandas as pd

import statsmodels.formula.api as smf
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

from pybedtools import BedTool
import matplotlib.pyplot as plt


def compute_gc_content(
    bin_info: pd.DataFrame, ref_file: str, mapp_file=None, genome_size=None
):
    """Compute per-bin GC fraction and optionally mean mappability.

    Uses pybedtools ``nucleotide_content`` on the reference FASTA. If a
    mappability BED track is provided, per-bin mean mappability is also computed.

    Parameters
    ----------
    bin_info : pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` columns defining genomic bins.
    ref_file : str
        Path to the reference genome FASTA.
    mapp_file : str or None
        Path to a BED-format mappability track (optional).
    genome_size : str or None
        Path to chromosome sizes file (required when *mapp_file* is given).

    Returns
    -------
    pd.DataFrame
        *bin_info* augmented with ``GC`` and ``MAP`` columns.
    """
    logging.info("compute GC content")
    gc_df = bin_info.merge(
        BedTool.from_dataframe(bin_info[["#CHR", "START", "END"]])
        .nucleotide_content(fi=ref_file)
        .to_dataframe(disable_auto_names=True)
        .rename(
            columns={
                "#1_usercol": "#CHR",
                "2_usercol": "START",
                "3_usercol": "END",
                "5_pct_gc": "GC",
            }
        )[["#CHR", "START", "END", "GC"]]
    )
    if mapp_file:
        logging.info("compute mappability")
        map_bt = BedTool(mapp_file)
        map_cov = (
            BedTool.from_dataframe(bin_info[["#CHR", "START", "END"]])
            .map(b=map_bt, c=4, o="mean", g=genome_size)
            .to_dataframe(disable_auto_names=True)
        )
        map_cov.columns = ["#CHR", "START", "END", "MAP"]
        map_cov["MAP"] = (
            pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(1.0).clip(0.0, 1.0)
        )
        gc_df = gc_df.merge(map_cov, on=["#CHR", "START", "END"], how="left")
    else:
        gc_df["MAP"] = 1.0
    return gc_df


def bias_correction_rdr(
    raw_rdr_mat: np.ndarray,
    gc_df: pd.DataFrame,
    rep_ids: list,
    has_mapp=False,
    out_dir=None,
    eps_quantile=0.01,
    gc_quantile=[0.01, 0.99],
):
    """Correct read-depth ratios (RDR) for GC-content (and optionally mappability) bias.

    Fits a median quantile regression of RDR on GC (quadratic) per tumor sample,
    divides raw RDR by the fitted expected RDR, and normalizes by the mean
    corrected RDR. Saves diagnostic scatter and hexbin plots.

    Parameters
    ----------
    raw_rdr_mat : np.ndarray
        Raw RDR matrix (bins x samples).
    gc_df : pd.DataFrame
        DataFrame with ``GC`` and ``MAP`` columns per bin.
    rep_ids : list[str]
        Replicate identifiers (tumor samples only).
    has_mapp : bool
        If True, include mappability as a covariate.
    out_dir : str or None
        Directory for diagnostic plots.
    eps_quantile : float
        Lower quantile of expected RDR used as floor to avoid division by zero.
    gc_quantile : list[float]
        Lower and upper GC quantiles defining the fitting range.

    Returns
    -------
    np.ndarray
        GC-corrected RDR matrix (bins x samples).
    """
    logging.info("correct for GC biases on RDR")
    gc = gc_df["GC"].to_numpy()
    mapv = gc_df["MAP"].to_numpy()

    gc_lo, gc_hi = np.nanquantile(gc, gc_quantile)
    logging.info(f"GC-content range for fitting: [{gc_lo}, {gc_hi}]")
    fit_mask = (gc >= gc_lo) & (gc <= gc_hi)

    gccorr_rdr_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        raw_rdrs = raw_rdr_mat[:, si]
        gc_df = pd.DataFrame({"RD": raw_rdrs[fit_mask], "GC": gc[fit_mask]})
        if has_mapp:
            gc_df["MAP"] = mapv[fit_mask]
            res = smf.quantreg("RD ~ GC + I(GC**2) + MAP + I(MAP**2)", data=gc_df).fit(
                q=0.5
            )
            logging.info(res.summary())
            exp_rdrs = res.predict(pd.DataFrame({"GC": gc, "MAP": mapv}))
        else:
            res = smf.quantreg("RD ~ GC + I(GC**2)", data=gc_df).fit(q=0.5)
            logging.info(res.summary())
            exp_rdrs = res.predict(pd.DataFrame({"GC": gc})).to_numpy()
        eps = np.nanquantile(exp_rdrs, eps_quantile)
        logging.info(f"inferred epsilon for exp RDR={eps}")
        den = exp_rdrs.copy()
        den[~np.isfinite(den)] = np.nan
        den = np.where(np.isfinite(den), den, eps)
        den = np.clip(den, eps, None)
        corr_rdrs = raw_rdrs / den

        # plot fitted line
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        x = np.linspace(gc.min(), gc.max(), 300)
        plt.scatter(gc, raw_rdrs, s=2, alpha=0.2)
        plt.plot(x, res.predict({"GC": x}), linewidth=2, color="red")
        plt.xlabel("GC")
        plt.ylabel("RDR")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{rep_id}.gc_corr_scatter.png"), dpi=300)
        plt.close(fig)

        corr_factor = np.mean(corr_rdrs)
        logging.info(f"GC correction factor={corr_factor}")

        # raw RDRs
        logging.info("hist: raw RDR")
        counts, edges = np.histogram(raw_rdrs)
        logging.info("bin_left\tbin_right\tcount")
        for l, r, c in zip(edges[:-1], edges[1:], counts):
            logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")

        # expected RDRs
        logging.info("hist: expected RDR corr_fit")
        counts, edges = np.histogram(exp_rdrs)
        logging.info("bin_left\tbin_right\tcount")
        for l, r, c in zip(edges[:-1], edges[1:], counts):
            logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")

        logging.info("hist: corr_rdrs")
        counts, edges = np.histogram(corr_rdrs)
        logging.info("bin_left\tbin_right\tcount")
        for l, r, c in zip(edges[:-1], edges[1:], counts):
            logging.info(f"{l:0.2f}\t{r:0.2f}\t{int(c)}")

        gccorr_rdr_mat[:, si] = corr_rdrs / corr_factor
        plot_gc_bias(
            gc,
            raw_rdr_mat[:, si],
            gccorr_rdr_mat[:, si],
            rep_id=rep_id,
            out_dir=out_dir,
        )
    return gccorr_rdr_mat


def plot_gc_bias(gc, raw_rdr, corr_rdr, rep_id=None, out_dir=None):
    """Produce a side-by-side hexbin plot comparing uncorrected vs. GC-corrected RDR.

    Each panel is annotated with the median absolute deviation (MAD).

    Parameters
    ----------
    gc : np.ndarray
        Per-bin GC content.
    raw_rdr : np.ndarray
        Uncorrected RDR values.
    corr_rdr : np.ndarray
        GC-corrected RDR values.
    rep_id : str or None
        Replicate identifier (used in the output filename).
    out_dir : str or None
        Directory for the output PNG.
    """
    def mad(x):
        return np.median(np.abs(x - np.median(x)))

    mad_raw = mad(raw_rdr)
    mad_corr = mad(corr_rdr)

    fig, axes = plt.subplots(1, 2, figsize=(6, 3), sharey=True)

    # --- Uncorrected ---
    hb1 = axes[0].hexbin(
        gc,
        raw_rdr,
        gridsize=100,
        cmap="Blues",
        mincnt=1,
        linewidths=0,
        reduce_C_function=np.median,
    )
    axes[0].set_xlabel("GC content")
    axes[0].set_ylabel("RDR")
    axes[0].set_title(f"Uncorrected RDR\nMAD={mad_raw:.3f}")

    # --- Corrected ---
    hb2 = axes[1].hexbin(
        gc,
        corr_rdr,
        gridsize=100,
        cmap="Blues",
        mincnt=1,
        linewidths=0,
        reduce_C_function=np.median,
    )
    axes[1].set_xlabel("GC content")
    axes[1].set_ylabel("RDR")
    axes[1].set_title(f"GC-corrected RDR\nMAD={mad_corr:.3f}")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{rep_id}.gc_corr.png"), dpi=300)
    plt.close(fig)
