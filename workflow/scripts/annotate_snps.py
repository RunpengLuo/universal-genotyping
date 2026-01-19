import os, sys, gzip, argparse
from pathlib import Path

import numpy as np
import pandas as pd

from snakemake.script import snakemake as sm

from utils import read_VCF

raw_snps = read_VCF(sm.input["raw_snp_file"])
chroms = sm.config["chromosomes"]

print("start annotate_snps")
assert all(c in raw_snps.columns for c in ["AD", "DP", "OTH"]), (
    "invalid cellsnp-lite format"
)
raw_snps["AD"] = raw_snps["AD"].astype(int)
raw_snps["DP"] = raw_snps["DP"].astype(int)
raw_snps["OTH"] = raw_snps["OTH"].astype(int)
raw_snps["REF_COUNT"] = raw_snps["DP"] - (raw_snps["AD"] + raw_snps["OTH"])

# filter sex-chromosome SNPs
raw_snps = raw_snps[
    raw_snps["#CHR"].astype(str).isin([f"chr{chrname}" for chrname in chroms])
]

raw_snps["VAF"] = raw_snps["AD"] / raw_snps["DP"]
raw_snps["REF_AGG"] = raw_snps["DP"] - raw_snps["AD"]


min_het_reads = sm.config["params_annotate_snps"]["min_het_reads"]
min_hom_dp = sm.config["params_annotate_snps"]["min_hom_dp"]
min_vaf_thres = sm.config["params_annotate_snps"]["min_vaf_thres"]
het = (
    (raw_snps["AD"] >= min_het_reads)
    & (raw_snps["REF_AGG"] >= min_het_reads)
    & (raw_snps["VAF"].between(min_vaf_thres, 1 - min_vaf_thres))
)
hom_alt = (raw_snps["AD"] == raw_snps["DP"]) & (raw_snps["DP"] >= min_hom_dp)
hom_ref = (raw_snps["AD"] == 0) & (raw_snps["DP"] >= min_hom_dp)

keep = het | hom_alt | hom_ref
print(
    f"#SNPs={keep.sum()}/{len(raw_snps)}, het={het.sum()}, hom_alt={hom_alt.sum()}, "
    f"hom_ref={hom_ref.sum()}"
)

raw_snps = raw_snps[keep].copy()

# GT
gt = np.full(len(raw_snps), "0/0", dtype=object)
gt[hom_alt.loc[raw_snps.index]] = "1/1"
gt[het.loc[raw_snps.index]] = "0/1"
raw_snps["FORMAT"] = "GT"
raw_snps["SAMPLE"] = gt

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
raw_snps = raw_snps[cols]

for chrname in chroms:
    chrom = f"chr{chrname}"
    out_snp_file = f"genotype/{chrom}.vcf"
    with open(out_snp_file, "w") as fd:
        fd.write("##fileformat=VCFv4.2\n")
        fd.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Pseudobulk genotype">\n'
        )
        fd.write("#" + "\t".join(cols) + "\n")
        raw_snps[raw_snps["CHROM"] == chrom].to_csv(
            fd, sep="\t", index=False, header=False
        )
print("finished annotate_snps")
