#!/usr/bin/env python

import os, gzip, argparse
from pathlib import Path

import numpy as np
import pandas as pd

from utils import read_barcodes, load_cellsnp_files


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--cellsnplite_result_dir", required=True)
    p.add_argument("-o", "--out_vcf", required=True)
    p.add_argument("-v", "--vaf_threshold", type=float, default=0.1)
    p.add_argument("--min_het_reads", type=int, default=2)
    p.add_argument("--min_hom_dp", type=int, default=10)
    a = p.parse_args()

    cellsnp_dir = os.path.abspath(a.cellsnplite_result_dir)

    # metadata
    base = list(Path(cellsnp_dir).glob("cellSNP.base*"))[0]
    meta = pd.read_csv(
        base,
        sep="\t",
        comment="#",
        usecols=["#CHR", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"],
    )
    print(f"#SNPs from cellsnp-lite genotyping={len(meta)}")

    # counts
    barcodes = read_barcodes(os.path.join(cellsnp_dir, "cellSNP.samples.tsv"))
    snp_df, dp_mat, ref_mat, alt_mat = load_cellsnp_files(cellsnp_dir, barcodes)

    snp_df = snp_df.merge(meta, on=["#CHR", "POS"], how="left")
    snp_df["DP"] = np.asarray(dp_mat.sum(axis=1)).ravel()
    snp_df["AD"] = np.asarray(alt_mat.sum(axis=1)).ravel()
    snp_df["REF_COUNT"] = np.asarray(ref_mat.sum(axis=1)).ravel()
    snp_df = snp_df[snp_df.DP > 0].copy()
    snp_df = snp_df[snp_df["#CHR"].astype(str).isin([str(i) for i in range(1, 23)])]
    snp_df["CHROM"] = "chr" + snp_df["#CHR"].astype(str)

    # filters
    snp_df["VAF"] = snp_df["AD"] / snp_df["DP"]
    snp_df["REF_AGG"] = snp_df["DP"] - snp_df["AD"]

    het = (
        (snp_df["AD"] >= a.min_het_reads)
        & (snp_df["REF_AGG"] >= a.min_het_reads)
        & (snp_df["VAF"].between(a.vaf_threshold, 1 - a.vaf_threshold))
    )
    hom_alt = (snp_df["AD"] == snp_df["DP"]) & (snp_df["DP"] >= a.min_hom_dp)
    hom_ref = (snp_df["AD"] == 0) & (snp_df["DP"] >= a.min_hom_dp)

    keep = het | hom_alt | hom_ref
    print(
        f"#SNPs={keep.sum()}/{len(snp_df)}, het={het.sum()}, hom_alt={hom_alt.sum()}, "
        f"hom_ref={hom_ref.sum()}"
    )

    snp_df = snp_df[keep].copy()

    # GT
    gt = np.full(len(snp_df), "0/0", dtype=object)
    gt[hom_alt.loc[snp_df.index]] = "1/1"
    gt[het.loc[snp_df.index]] = "0/1"
    snp_df["SAMPLE"] = gt

    # collapse duplicates
    snp_df["KEY"] = (
        snp_df["CHROM"]
        + "_"
        + snp_df["POS"].astype(str)
        + "_"
        + snp_df["REF"].astype(str)
        + "_"
        + snp_df["ALT"].astype(str)
    )

    g = (
        snp_df.groupby("KEY", as_index=False)
        .agg(
            CHROM=("CHROM", "first"),
            POS=("POS", "first"),
            ID=("ID", "first"),
            REF=("REF", "first"),
            ALT=("ALT", "first"),
            QUAL=("QUAL", "first"),
            FILTER=("FILTER", "first"),
            AD=("AD", "sum"),
            DP=("DP", "sum"),
            REF_COUNT=("REF_COUNT", "sum"),
            SAMPLE=("SAMPLE", lambda x: x.value_counts().index[0]),
        )
        .sort_values(["CHROM", "POS"])
    )

    g["OTH"] = (g["DP"] - g["AD"] - g["REF_COUNT"]).clip(lower=0)
    g["INFO"] = (
        "AD="
        + g["AD"].astype(str)
        + ";DP="
        + g["DP"].astype(str)
        + ";OTH="
        + g["OTH"].astype(str)
    )
    g["FORMAT"] = "GT"

    cols = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SAMPLE",
    ]
    g = g[cols]

    out_dir = os.path.dirname(a.out_vcf)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with gzip.open(a.out_vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Pseudobulk genotype">\n'
        )
        fh.write("#" + "\t".join(cols) + "\n")
        g.to_csv(fh, sep="\t", index=False, header=False)
