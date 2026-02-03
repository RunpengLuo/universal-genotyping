import os, sys, gzip, logging, subprocess
from snakemake.script import snakemake as sm

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import scipy
import numpy as np
import pandas as pd

from utils import read_VCF

##################################################
"""
Annotate genotpyed SNPs from pseudobulk samples.
Merge SNPs across VCFs from multiple sources.

After annotation, only unique bi-allelic Het or Hom-Alt SNPs are kept.
"""
logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logging.info("start annotate_snps")

chroms = sm.config["chromosomes"]
filter_nz_OTH = sm.params["filter_nz_OTH"]
filter_hom_ALT = sm.params["filter_hom_ALT"]
min_het_reads = sm.params["min_het_reads"]
min_hom_dp = sm.params["min_hom_dp"]
min_vaf_thres = sm.params["min_vaf_thres"]
snp_lists = []
KEY = ["#CHROM", "POS", "REF", "ALT"]
CNT = ["DP", "AD", "OTH"]


def get_gt(row):
    if row["is_het"]:
        return "0/1"
    elif row["is_hom_alt"]:
        return "1/1"
    elif row["is_hom_ref"]:
        return "0/0"
    else:
        return "./."


keep_gts = ["0/1"]
if not filter_hom_ALT:
    keep_gts.append("1/1")

##################################################
for data_type, raw_snp_file in zip(sm.params["data_types"], sm.input["raw_snp_files"]):
    raw_snps = read_VCF(raw_snp_file)
    assert all(c in raw_snps.columns for c in CNT), "invalid cellsnp-lite format"
    num_raw_snps = len(raw_snps)
    raw_snps["SOURCE"] = data_type
    raw_snps["FORMAT"] = "GT"
    raw_snps["AD"] = raw_snps["AD"].astype(int)
    raw_snps["DP"] = raw_snps["DP"].astype(int)
    raw_snps["OTH"] = raw_snps["OTH"].astype(int)
    raw_snps["REF_COUNT"] = raw_snps["DP"] - raw_snps["AD"]
    raw_snps["AR"] = raw_snps["AD"] / raw_snps["DP"]

    # genotyping criteria
    raw_snps["is_het"] = (
        raw_snps[["AD", "REF_COUNT"]].min(axis=1) >= min_het_reads
    ) & raw_snps["AR"].between(min_vaf_thres, 1 - min_vaf_thres)
    raw_snps["is_hom_alt"] = (raw_snps["AD"] == raw_snps["DP"]) & (
        raw_snps["DP"] >= min_hom_dp
    )
    raw_snps["is_hom_ref"] = (raw_snps["AD"] == 0) & (raw_snps["DP"] >= min_hom_dp)
    raw_snps["SAMPLE"] = raw_snps.apply(get_gt, axis=1)

    raw_snps = raw_snps[raw_snps["SAMPLE"].isin(keep_gts)].reset_index(drop=True)

    # filter SNPs by chromosomes
    raw_snps = raw_snps[
        raw_snps["#CHROM"].astype(str).isin([f"chr{chrname}" for chrname in chroms])
    ]

    # filter multi-allelic SNPs, ~0.1%
    if filter_nz_OTH:
        nz_oth = raw_snps["OTH"] > 0
        n_nz = int(nz_oth.sum())

        n_filt_hom_alt = int(raw_snps.loc[nz_oth, "is_hom_alt"].sum())
        n_filt_het = int(raw_snps.loc[nz_oth, "is_het"].sum())
        n_filt_hom_ref = int(raw_snps.loc[nz_oth, "is_hom_ref"].sum())
        n_filt_other = n_nz - (n_filt_hom_alt + n_filt_het + n_filt_hom_ref)

        logging.info(
            f"{data_type}, #nz-OTH SNPs={n_nz}/{num_raw_snps} "
            f"(hom_alt={n_filt_hom_alt}, het={n_filt_het}, hom_ref={n_filt_hom_ref}, other={n_filt_other})"
        )

        raw_snps = raw_snps[~nz_oth]

    # filter duplicated SNPs
    raw_snps = raw_snps.drop_duplicates(subset=KEY, keep="first").reset_index(drop=True)

    num_het = (raw_snps["SAMPLE"] == "0/1").sum()
    num_hom_ref = (raw_snps["SAMPLE"] == "0/0").sum()
    num_hom_alt = (raw_snps["SAMPLE"] == "1/1").sum()
    logging.info(
        f"{data_type}, #kept SNPs={len(raw_snps)}/{num_raw_snps}, #het={num_het}, #hom_alt={num_hom_alt}, #hom_ref={num_hom_ref}"
    )
    snp_lists.append(raw_snps)

##################################################
final_snps = snp_lists[0]
if len(snp_lists) > 1:
    # merge SNP file over multiple sources
    final_snps = pd.concat(snp_lists, axis=0, ignore_index=True)

    # 1. filter duplicated SNPs common in multiple sources
    final_snps = final_snps.drop_duplicates(subset=KEY, keep="first").reset_index(
        drop=True
    )

    # 2. filter multi-allelic SNPs.
    is_dup = final_snps.duplicated(subset=["#CHROM", "POS"], keep=False)
    final_snps = final_snps[~is_dup].reset_index(drop=True)

final_snps["#CHROM"] = pd.Categorical(
    final_snps["#CHROM"].astype(str),
    categories=[f"chr{c}" for c in chroms],
    ordered=True,
)
final_snps = final_snps.sort_values(["#CHROM", "POS"], kind="mergesort").reset_index(
    drop=True
)

cols = [
    "#CHROM",
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
final_snps = final_snps[cols]

final_snps_chs = final_snps.groupby("#CHROM", sort=False)
for chrname, out_snp_file in zip(chroms, sm.output["snp_files"]):
    chrom = f"chr{chrname}"
    with open(out_snp_file[:-3], "w") as fd:
        fd.write("##fileformat=VCFv4.2\n")
        fd.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Pseudobulk genotype">\n'
        )
        fd.write("\t".join(cols) + "\n")
        if chrom in final_snps_chs.groups:
            final_snps_chs.get_group(chrom)[cols].to_csv(
                fd, sep="\t", index=False, header=False
            )

    # out_snp_file[:-3] will be removed by bgzip.
    subprocess.run(["bgzip", "-f", out_snp_file[:-3]], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", out_snp_file], check=True)
logging.info("finished annotate_snps")
