import os
import sys
import subprocess
from io import StringIO

import pandas as pd
import numpy as np


def read_VCF(vcf_file: str, addchr=True, addkey=False):
    """
    load VCF file as dataframe.
    If phased, parse GT[0] as USEREF, check PS
    """
    snps = pd.read_csv(vcf_file, comment="#", sep="\t", header=None)
    if snps.empty:
        return None
    ncols = snps.shape[1]
    assert ncols == 8 or ncols >= 10, "invalid VCF file"
    colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if ncols >= 10:
        colnames += ["FORMAT", "SAMPLE"]
        snps = snps.iloc[:, :10].copy()
    snps.columns = colnames

    if addchr and not str(snps["#CHROM"].iloc[0]).startswith("chr"):
        snps["#CHROM"] = "chr" + snps["#CHROM"].astype(str)
    snps["#CHR"] = snps["#CHROM"]

    # parse INFO column
    info_kvs = (
        snps["INFO"]
        .fillna("")
        .str.split(";")
        .explode()
        .loc[lambda s: s.ne("")]
        .to_frame("kv")
    )
    kv = info_kvs["kv"].str.split("=", n=1, expand=True)
    info_kvs["key"] = kv[0]
    info_kvs["val"] = kv[1] if kv.shape[1] > 1 else None
    info_kvs["val"] = info_kvs["val"].fillna(True)  # INFO flags (no '=') -> True
    info_kvs["row"] = info_kvs.index

    info_wide = info_kvs.pivot_table(
        index="row", columns="key", values="val", aggfunc="first"
    )
    snps = snps.join(info_wide)

    # parse FORMAT column
    if "FORMAT" in snps.columns:
        fmt_keys = snps["FORMAT"].fillna("").str.split(":")
        samp_vals = snps["SAMPLE"].fillna("").str.split(":")

        fmt_long = pd.DataFrame(
            {"key": fmt_keys.explode(), "val": samp_vals.explode()}
        ).dropna(subset=["key"])
        fmt_long["row"] = fmt_keys.explode().index  # original variant row index

        fmt_wide = fmt_long.pivot_table(
            index="row", columns="key", values="val", aggfunc="first"
        )
        snps = snps.drop(columns=["FORMAT", "SAMPLE"]).join(fmt_wide)
    
    if "GT" in snps.columns and "PS" not in snps.columns:
        # global phased block if PS information N/A
        snps["PS"] = 1
    if addkey:
        snps["KEY"] = snps["#CHROM"].astype(str) + "_" + snps["POS"].astype(str)
    snps = snps.reset_index(drop=True)
    return snps

def symlink_force(src, dst):
    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    os.symlink(os.path.abspath(src), os.path.abspath(dst))
