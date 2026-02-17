import os
import sys
import subprocess
from io import StringIO
from collections import OrderedDict

import pandas as pd
import numpy as np
import pyranges as pr

from utils import *


def get_chr_sizes(sz_file: str):
    chr_sizes = OrderedDict()
    with open(sz_file, "r") as rfd:
        for line in rfd.readlines():
            ch, sizes = line.strip().split()
            chr_sizes[ch] = int(sizes)
        rfd.close()
    return chr_sizes


def read_VCF(
    vcf_file: str,
    addchr=True,
    addkey=False,
    snps_presorted=False,
    add_pos0=False,
    add_phase1=False,
):
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
    snps["POS"] = snps["POS"].astype(np.int64)
    snps["RAW_SNP_IDX"] = np.arange(len(snps))
    if addchr and not str(snps["#CHROM"].iloc[0]).startswith("chr"):
        snps["#CHROM"] = "chr" + snps["#CHROM"].astype(str)

    # sort snps by #CHROM and POS
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM", "chrMT"]
    snps["#CHROM"] = snps["#CHROM"].str.replace("^chrMT$", "chrM", regex=True)
    snps["#CHROM"] = pd.Categorical(
        snps["#CHROM"], categories=chrom_order, ordered=True
    )
    if not snps_presorted:
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

    if addkey:
        snps["KEY"] = snps["#CHROM"].astype(str) + "_" + snps["POS"].astype(str)
    if add_pos0:
        snps["POS0"] = snps["POS"] - 1
    if add_phase1:
        snps["PHASE"] = snps["GT"].str[2].astype(np.int8).to_numpy()
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


def read_ucn_file(seg_ucn_file: str):
    segs_df = pd.read_table(seg_ucn_file, sep="\t")
    segs_df = sort_df_chr(segs_df, pos="START")

    n_clones = len([cname for cname in segs_df.columns if cname.startswith("cn_")])
    clones = ["normal"] + [f"clone{c}" for c in range(1, n_clones)]
    segs_df["CNP"] = segs_df.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    segs_df["PROPS"] = segs_df.apply(
        func=lambda r: ";".join(str(r[f"u_{c}"]) for c in clones), axis=1
    )
    return segs_df, clones


def read_barcodes(bc_file: str):
    barcodes = (
        pd.read_table(bc_file, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()
    )
    return barcodes


def read_celltypes(celltype_file: str, ref_label):
    celltypes = pd.read_table(celltype_file, sep="\t")[["BARCODE", ref_label]].reset_index(drop=True)
    return celltypes


def read_genes_gtf_file(gtf_file: str, id_col="gene_ids"):
    gr = pr.read_gtf(gtf_file, rename_attr=True)
    genes = (
        gr.df.query("Feature == 'gene'")[["Chromosome", "Start", "End", "gene_id"]]
        .drop_duplicates("gene_id", keep="first")
        .rename(
            columns={
                "gene_id": id_col,
                "Chromosome": "#CHR",
                "Start": "START1",
                "End": "END",
            }
        )
    )
    # convert to 0-based BED format
    genes["START"] = genes["START1"] - 1

    return genes


def read_genes_bed_file(bed_file: str, id_col="gene_id"):
    gr = pr.read_bed(bed_file, as_df=True)
    genes = (
        gr[["Chromosome", "Start", "End", "Name"]]
        .drop_duplicates("Name", keep="first")
        .rename(
            columns={
                "Name": id_col,
                "Chromosome": "#CHR",
                "Start": "START",
                "End": "END",
            }
        )
    )
    return genes
