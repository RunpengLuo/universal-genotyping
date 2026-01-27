import os
import sys
import subprocess
from io import StringIO
from collections import OrderedDict

import pandas as pd
import numpy as np


def get_chr2ord(ch):
    chr2ord = {}
    for i in range(1, 23):
        chr2ord[f"{ch}{i}"] = i
    chr2ord[f"{ch}X"] = 23
    chr2ord[f"{ch}Y"] = 24
    chr2ord[f"{ch}M"] = 25
    return chr2ord


def sort_chroms(chromosomes: list):
    assert len(chromosomes) != 0
    chromosomes = [str(c) for c in chromosomes]
    ch = "chr" if str(chromosomes[0]).startswith("chr") else ""
    chr2ord = get_chr2ord(ch)
    return sorted(chromosomes, key=lambda x: chr2ord[x])


def sort_df_chr(df: pd.DataFrame, ch="#CHR", pos="POS"):
    chs = sort_chroms(df[ch].unique().tolist())
    df[ch] = pd.Categorical(df[ch], categories=chs, ordered=True)
    df.sort_values(by=[ch, pos], inplace=True, ignore_index=True)
    return df


def get_chr_sizes(sz_file: str):
    chr_sizes = OrderedDict()
    with open(sz_file, "r") as rfd:
        for line in rfd.readlines():
            ch, sizes = line.strip().split()
            chr_sizes[ch] = int(sizes)
        rfd.close()
    return chr_sizes


def read_VCF(vcf_file: str, addchr=True, addkey=False):
    """
    load VCF file as dataframe.
    If phased, parse GT[0] as USEREF, check PS
    """
    snps = pd.read_csv(
        vcf_file, comment="#", sep="\t", header=None, dtype={0: "string"}
    )
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

    # sort snps by #CHROM and POS
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM", "chrMT"]
    snps["#CHROM"] = snps["#CHROM"].str.replace("^chrMT$", "chrM", regex=True)
    snps["#CHROM"] = pd.Categorical(
        snps["#CHROM"], categories=chrom_order, ordered=True
    )
    snps = snps.sort_values(["#CHROM", "POS"], kind="mergesort")

    # copy #CHR field
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


def read_region_file(region_bed_file: str, addchr=True):
    regions = pd.read_table(
        region_bed_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["Chromosome", "Start", "End"],
        dtype={0: "string"},
    )
    if not str(regions["Chromosome"].iloc[0]).startswith("chr") and addchr:
        regions["Chromosome"] = "chr" + regions["Chromosome"].astype(str)
    regions["#CHR"] = regions["Chromosome"]
    regions["START"] = regions["Start"]
    regions["END"] = regions["End"]
    return regions


def symlink_force(src, dst):
    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    os.symlink(os.path.abspath(src), os.path.abspath(dst))
