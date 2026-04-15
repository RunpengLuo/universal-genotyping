import os
from collections import OrderedDict

import pandas as pd
import numpy as np
import pyranges as pr

from utils import *


def get_chr_sizes(sz_file: str):
    """Read a two-column chromosome-sizes file and return an OrderedDict mapping name to length.

    Parameters
    ----------
    sz_file : str
        Path to a tab-separated file with columns (chromosome, size).

    Returns
    -------
    OrderedDict[str, int]
        Chromosome name to integer length.
    """
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

    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM", "chrMT"]
    snps["#CHROM"] = snps["#CHROM"].str.replace("^chrMT$", "chrM", regex=True)
    snps["#CHROM"] = pd.Categorical(
        snps["#CHROM"], categories=chrom_order, ordered=True
    )
    if not snps_presorted:
        snps = snps.sort_values(["#CHROM", "POS"], kind="mergesort")

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
    """Read the first three columns of a BED file into a DataFrame.

    Parameters
    ----------
    region_bed_file : str
        Path to a BED file.
    addchr : bool
        If True, prefix chromosome names with ``chr`` when missing.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` (and pyranges-compatible aliases).
    """
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
    """Read a HATCHet-style UCN segment file and derive CNP and PROPS columns.

    Parameters
    ----------
    seg_ucn_file : str
        Path to a tab-separated UCN file with per-clone copy-number columns.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        Sorted segment DataFrame with ``CNP`` and ``PROPS`` columns, and
        the list of clone names.
    """
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
    """Read a barcode file (one barcode per line) and return as a list of strings.

    Parameters
    ----------
    bc_file : str
        Path to a text file with one barcode per line.

    Returns
    -------
    list[str]
        Barcodes.
    """
    barcodes = (
        pd.read_table(bc_file, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()
    )
    return barcodes


def compute_depth_statistics(dp_raw, win_df, sample_ids, rep_ids):
    """Compute per-chromosome and whole-genome mean/median depth per sample.

    Returns a DataFrame with columns: SAMPLE, #CHR, REP_ID, mean_depth, median_depth.
    """
    chroms = win_df["#CHR"].to_numpy()
    chrom_order = {c: i for i, c in enumerate(CHROM_ORDER)}
    sorted_chroms = sorted(
        win_df["#CHR"].unique(),
        key=lambda c: chrom_order.get(c, len(chrom_order)),
    )
    rows = []
    for chrom in sorted_chroms:
        mask = chroms == chrom
        for s, rep_id in enumerate(rep_ids):
            vals = dp_raw[mask, s]
            rows.append([sample_ids[s], chrom, rep_id, float(np.mean(vals)), float(np.median(vals))])
    for s, rep_id in enumerate(rep_ids):
        vals = dp_raw[:, s]
        rows.append([sample_ids[s], "TOTAL", rep_id, float(np.mean(vals)), float(np.median(vals))])
    return pd.DataFrame(rows, columns=["SAMPLE", "#CHR", "REP_ID", "mean_depth", "median_depth"])


def read_genes_gtf_file(gtf_file: str, id_col="gene_ids"):
    """Parse a GTF file and return gene-level records with genomic coordinates.

    Extracts gene features, deduplicates by ``gene_id``, and converts to
    0-based BED-like coordinates.

    Parameters
    ----------
    gtf_file : str
        Path to a GTF annotation file.
    id_col : str
        Column name for the gene identifier in the output.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``#CHR``, ``START`` (0-based), ``END``, and *id_col*.
    """
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


def read_exons_gtf_file(gtf_file: str):
    """Parse a GTF file and return exon-level records with genomic coordinates.

    Extracts exon features and converts to 0-based BED-like coordinates.

    Parameters
    ----------
    gtf_file : str
        Path to a GTF annotation file.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``#CHR``, ``START`` (0-based), ``END``, and ``gene_id``.
    """
    gr = pr.read_gtf(gtf_file, rename_attr=True)
    exons = gr.df.query("Feature == 'exon'")[
        ["Chromosome", "Start", "End", "gene_id"]
    ].rename(
        columns={
            "Chromosome": "#CHR",
            "Start": "START1",
            "End": "END",
        }
    )
    exons["START"] = exons["START1"] - 1
    return exons[["#CHR", "START", "END", "gene_id"]]


def read_genes_bed_file(bed_file: str, id_col="gene_id"):
    """Read a BED file and return deduplicated gene-level records.

    Parameters
    ----------
    bed_file : str
        Path to a BED file with gene names in the Name column.
    id_col : str
        Column name for the gene identifier in the output.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END``, and *id_col*.
    """
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
