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

from io_utils import read_VCF

##################################################
"""
Given cellsnp-lite results from multiple replicates.
After annotation, only unique bi-allelic Het or Hom-Alt SNPs are kept.
"""
logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

chroms = sm.config["chromosomes"]
filter_nz_OTH = sm.params["filter_nz_OTH"]
filter_hom_ALT = sm.params["filter_hom_ALT"]
min_het_reads = sm.params["min_het_reads"]
min_hom_dp = sm.params["min_hom_dp"]
min_vaf_thres = sm.params["min_vaf_thres"]
snp_lists = []
KEY = ["#CHROM", "POS", "REF", "ALT"]
CNT = ["DP", "AD", "OTH"]

logging.info(
    f"start annotate_snps, filter_nz_OTH={filter_nz_OTH}, filter_hom_ALT={filter_hom_ALT}"
)
logging.info(
    f"min_het_reads={min_het_reads}, min_hom_dp={min_hom_dp}, min_vaf_thres={min_vaf_thres}"
)


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
data_types = list(sm.params["data_types"])
rep_ids = list(sm.params["rep_ids"])
raw_snp_files = list(sm.params["raw_snp_files"])
assert len(data_types) == 1, "multimodal genotyping unsupported yet"
data_type = data_types[0]

# step 1. remove duplicated rows and collapse VCF files from multiple replicates
raw_snps_list = []
for idx, rep_id in enumerate(rep_ids):
    raw_snps = read_VCF(raw_snp_files[idx])
    raw_snps["KEY"] = (
        raw_snps["#CHROM"].astype(str)
        + "_"
        + raw_snps["POS"].astype(str)
        + "_"
        + raw_snps["REF"].astype(str)
        + "_"
        + raw_snps["ALT"].astype(str)
    )
    assert all(c in raw_snps.columns for c in CNT), "invalid cellsnp-lite format"
    for cnt in CNT:
        raw_snps[f"{cnt}{idx}"] = raw_snps[cnt].astype(np.int64)

    # filter SNPs by chrnames and duplicated entries.
    raw_snps = raw_snps[
        raw_snps["#CHROM"].astype(str).isin([f"chr{chrname}" for chrname in chroms])
    ]

    # filter rows with duplicated snp ids
    n_dup_rows = raw_snps.duplicated(subset="KEY", keep=False).sum()
    n_dup_keys = raw_snps.loc[
        raw_snps.duplicated(subset="KEY", keep=False), "KEY"
    ].nunique()
    logging.info(
        f"[{data_type} {rep_id}] duplicated rows, key=#CHROM_POS_REF_ALT: rows={n_dup_rows}, keys={n_dup_keys}"
    )
    raw_snps = raw_snps.drop_duplicates(subset="KEY", keep="first").reset_index(
        drop=True
    )
    raw_snps_list.append(raw_snps)

# union over all replicates
base_snps = pd.concat(
    [raw_snps[["KEY", "#CHROM", "POS", "REF", "ALT"]] for raw_snps in raw_snps_list],
    ignore_index=True,
).drop_duplicates(subset="KEY", keep="first")
base_snps = base_snps.sort_values(["#CHROM", "POS"], kind="mergesort")
for cnt in CNT:
    base_snps[cnt] = 0
for idx, raw_snps in enumerate(raw_snps_list):
    right_cols = ["KEY"] + [f"{cnt}{idx}" for cnt in CNT]
    base_snps = pd.merge(
        left=base_snps,
        right=raw_snps[right_cols],
        on="KEY",
        how="left",
        validate="one_to_one",
        sort=False,
    )
    for cnt in CNT:
        base_snps[cnt] += base_snps[f"{cnt}{idx}"].fillna(0).astype(np.int64)
    base_snps = base_snps.drop(columns=[f"{cnt}{idx}" for cnt in CNT])

# filter multi-allelic SNPs.
dup_mask = base_snps.duplicated(subset=["#CHROM", "POS"], keep=False)
n_dup_rows = int(dup_mask.sum())
n_dup_keys = int(base_snps.loc[dup_mask, ["#CHROM", "POS"]].drop_duplicates().shape[0])
logging.info(
    f"[{data_type}] duplicated rows, key=#CHROM_POS: rows={n_dup_rows}, keys={n_dup_keys}"
)
base_snps = base_snps.loc[~dup_mask, :].reset_index(drop=True)
nsnps_before_genotyping = len(base_snps)

# SNP genotyping.
base_snps["REF_COUNT"] = base_snps["DP"] - base_snps["AD"]
base_snps["AR"] = base_snps["AD"] / base_snps["DP"]

# genotyping criteria
base_snps["is_het"] = (
    base_snps[["AD", "REF_COUNT"]].min(axis=1) >= min_het_reads
) & base_snps["AR"].between(min_vaf_thres, 1 - min_vaf_thres)
base_snps["is_hom_alt"] = (base_snps["AD"] == base_snps["DP"]) & (
    base_snps["DP"] >= min_hom_dp
)
base_snps["is_hom_ref"] = (base_snps["AD"] == 0) & (base_snps["DP"] >= min_hom_dp)
base_snps["SAMPLE"] = base_snps.apply(get_gt, axis=1)

# filter SNPs by genotypes
base_snps = base_snps[base_snps["SAMPLE"].isin(keep_gts)].reset_index(drop=True)

nz_oth = base_snps["OTH"] > 0
n_nz = np.sum(nz_oth)

nz_oth_hom_alt = int(base_snps.loc[nz_oth, "is_hom_alt"].sum())
nz_oth_het = int(base_snps.loc[nz_oth, "is_het"].sum())
nz_oth_hom_ref = int(base_snps.loc[nz_oth, "is_hom_ref"].sum())
nz_oth_other = n_nz - (nz_oth_hom_alt + nz_oth_het + nz_oth_hom_ref)

logging.info(
    f"{data_type}, #nz-OTH SNPs={n_nz}/{nsnps_before_genotyping} "
    f"(hom_alt={nz_oth_hom_alt}, het={nz_oth_het}, hom_ref={nz_oth_hom_ref}, other={nz_oth_other})"
)
if filter_nz_OTH:
    base_snps = base_snps[~nz_oth]

num_het = (base_snps["SAMPLE"] == "0/1").sum()
num_hom_ref = (base_snps["SAMPLE"] == "0/0").sum()
num_hom_alt = (base_snps["SAMPLE"] == "1/1").sum()
num_snp_kept = len(base_snps)
logging.info(
    f"{data_type}, #kept SNPs={num_snp_kept}/{nsnps_before_genotyping}, #het={num_het}, #hom_alt={num_hom_alt}, #hom_ref={num_hom_ref}"
)

final_snps = base_snps.sort_values(["#CHROM", "POS"], kind="mergesort").reset_index(
    drop=True
)

final_snps["FORMAT"] = "GT"
final_snps["ID"] = "."
final_snps["QUAL"] = "."
final_snps["FILTER"] = "PASS"
final_snps["INFO"] = (
    "AD="
    + final_snps["AD"].astype(str)
    + ";DP="
    + final_snps["DP"].astype(str)
    + ";OTH="
    + final_snps["OTH"].astype(str)
)

##################################################
# save final vcf file
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
        fd.writelines(
            [
                "##fileformat=VCFv4.2\n",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Pseudobulk genotype">\n',
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth, REF+ALT">\n',
                '##INFO=<ID=AD,Number=1,Type=Integer,Description="Allele Depth for ALT allele">\n',
                '##INFO=<ID=OTH,Number=1,Type=Integer,Description="Allele Depth other than REF and ALT">\n',
            ]
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
