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
    gc_df["GC"] /= 100  # convert to %GC
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
    raw_rdr_mat: np.ndarray, gc_df: pd.DataFrame, rep_ids: list, out_dir=None
):
    logging.info("correct for GC biases on RDR")
    gc = gc_df["GC"].to_numpy()
    mapv = gc_df["MAP"].to_numpy()
    gccorr_rdr_mat = np.zeros_like(raw_rdr_mat, dtype=np.float32)
    for si, rep_id in enumerate(rep_ids):
        gc_df = pd.DataFrame({"RD": raw_rdr_mat[:, si], "GC": gc, "MAP": mapv})
        if np.any(gc_df["MAP"] != 1):
            # mappability
            mod = smf.quantreg("RD ~ GC + I(GC**2) + MAP + I(MAP**2)", data=gc_df).fit(
                q=0.5
            )
            gc_df["CORR_FIT"] = mod.predict(gc_df[["GC", "MAP"]])
            corr_rdrs = gc_df["RD"] / gc_df["CORR_FIT"].where(
                (gc_df["CORR_FIT"] > 0) & ~pd.isnull(gc_df["CORR_FIT"]), 1
            )
        else:
            mod = smf.quantreg("RD ~ GC + I(GC**2)", data=gc_df).fit(q=0.5)
            gc_df["CORR_FIT"] = mod.predict(gc_df[["GC"]])
            corr_rdrs = gc_df["RD"] / gc_df["CORR_FIT"].where(
                (gc_df["CORR_FIT"] > 0) & ~pd.isnull(gc_df["CORR_FIT"]), 1
            )
        corr_factor = np.mean(corr_rdrs)
        logging.info(f"GC correction factor={corr_factor}")
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
