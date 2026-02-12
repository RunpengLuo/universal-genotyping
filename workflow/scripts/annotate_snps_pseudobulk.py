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

from io_utils import read_VCF, get_chr_sizes

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
min_het_reads = int(sm.params["min_het_reads"])
min_hom_dp = int(sm.params["min_hom_dp"])
min_vaf_thres = float(sm.params["min_vaf_thres"])

logging.info(
    f"start annotate_snps_pseudobulk, filter_nz_OTH={filter_nz_OTH}, filter_hom_ALT={filter_hom_ALT}"
)
logging.info(
    f"min_het_reads={min_het_reads}, min_hom_dp={min_hom_dp}, min_vaf_thres={min_vaf_thres}"
)


def get_genotype(row):
    if row["is_het"]:
        return "0/1"
    if row["is_hom_alt"]:
        return "1/1"
    if row["is_hom_ref"]:
        return "0/0"
    return "./."


##################################################
# step 1. remove duplicated rows and collapse VCF files from multiple modalities
KEY = ["#CHROM", "POS", "REF", "ALT"]
CNT = ["DP", "AD", "OTH"]
modalities = list(sm.params["modalities"])
raw_snp_vcfs = list(sm.input["raw_snp_vcfs"])
raw_snps_list = []
for idx, modality in enumerate(modalities):
    raw_snps = read_VCF(raw_snp_vcfs[idx], addkey=True)
    assert all(c in raw_snps.columns for c in KEY + CNT), "invalid cellsnp-lite format"
    for cnt in CNT:
        raw_snps[f"{cnt}{idx}"] = raw_snps[cnt].astype(np.int64)

    # filter SNPs by chrnames and duplicated entries.
    raw_snps = raw_snps[
        raw_snps["#CHROM"].astype(str).isin([f"chr{chrname}" for chrname in chroms])
    ]

    # filter rows with duplicated snp ids
    n_dup_rows = raw_snps.duplicated(subset="KEY", keep=False).sum()
    if n_dup_rows > 0:
        n_dup_keys = raw_snps.loc[
            raw_snps.duplicated(subset="KEY", keep=False), "KEY"
        ].nunique()
        logging.warning(
            f"{modality} have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
        )
        logging.warning("drop duplicated rows.")
        raw_snps = raw_snps.drop_duplicates(subset="KEY", keep="first").reset_index(
            drop=True
        )
    raw_snps_list.append(raw_snps)

## union over all replicates
base_snps = pd.concat(
    [raw_snps[["KEY", "#CHROM", "POS", "REF", "ALT"]] for raw_snps in raw_snps_list],
    ignore_index=True,
).drop_duplicates(subset=["#CHROM", "POS", "REF", "ALT"], keep="first")

## filter any SNP rows that present in multiple modalities but has different allele pair.
dup_mask = base_snps.duplicated(subset="KEY", keep=False)
n_dup_rows = int(dup_mask.sum())
if n_dup_rows > 0:
    n_dup_keys = int(base_snps.loc[dup_mask, "KEY"].drop_duplicates().shape[0])
    logging.warning(
        f"merged df have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
    )
    logging.warning("drop duplicated rows.")
    base_snps = base_snps.loc[~dup_mask, :]
base_snps = base_snps.reset_index(drop=True)
nsnps_before_genotyping = len(base_snps)

base_snps = base_snps.sort_values(["#CHROM", "POS"], kind="mergesort")
## aggregate DP/AD/OTH counts across modalities
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

##################################################
# step 2. genotype and filter SNPs
base_snps["REF_COUNT"] = base_snps["DP"] - base_snps["AD"]
base_snps["AR"] = base_snps["AD"] / base_snps["DP"]

## genotyping criteria
base_snps["is_het"] = (
    base_snps[["AD", "REF_COUNT"]].min(axis=1) >= min_het_reads
) & base_snps["AR"].between(min_vaf_thres, 1 - min_vaf_thres)

base_snps["is_hom_alt"] = (base_snps["AD"] == base_snps["DP"]) & (
    base_snps["DP"] >= min_hom_dp
)
base_snps["is_hom_ref"] = (base_snps["AD"] == 0) & (base_snps["DP"] >= min_hom_dp)

## assign&filter genotype
base_snps["SAMPLE"] = base_snps.apply(get_genotype, axis=1)
keep_gts = ["0/1"]
if not filter_hom_ALT:
    keep_gts.append("1/1")
base_snps = base_snps[base_snps["SAMPLE"].isin(keep_gts)].reset_index(drop=True)

## handle OTH rows
nz_oth = base_snps["OTH"] > 0
nnz_oth = np.sum(nz_oth)
nz_oth_hom_alt = int(base_snps.loc[nz_oth, "is_hom_alt"].sum())
nz_oth_het = int(base_snps.loc[nz_oth, "is_het"].sum())
nz_oth_hom_ref = int(base_snps.loc[nz_oth, "is_hom_ref"].sum())
nz_oth_other = nnz_oth - (nz_oth_hom_alt + nz_oth_het + nz_oth_hom_ref)
logging.info(
    f"#nz-OTH SNPs={nnz_oth}/{nsnps_before_genotyping} "
    f"(hom_alt={nz_oth_hom_alt}, het={nz_oth_het}, hom_ref={nz_oth_hom_ref}, other={nz_oth_other})"
)
if filter_nz_OTH:
    logging.info("SNPs with nonzero OTHs are filtered.")
    base_snps = base_snps[~nz_oth]

## final summary stats
num_het = (base_snps["SAMPLE"] == "0/1").sum()
num_hom_ref = (base_snps["SAMPLE"] == "0/0").sum()
num_hom_alt = (base_snps["SAMPLE"] == "1/1").sum()
num_snp_kept = len(base_snps)
logging.info(f"#kept SNPs={num_snp_kept}/{nsnps_before_genotyping}")
logging.info(f"#het={num_het}, #hom_alt={num_hom_alt}, #hom_ref={num_hom_ref}")

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
chrom_sizes = get_chr_sizes(sm.input["genome_size"])
for chrname, out_snp_vcf in zip(chroms, sm.output["snp_vcfs"]):
    chrom = f"chr{chrname}"
    chrom_length = chrom_sizes[chrom]
    with open(out_snp_vcf[:-3], "w") as fd:
        fd.writelines(
            [
                "##fileformat=VCFv4.2\n",
                f"##contig=<ID={chrom},length={chrom_length}>\n",
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

    # out_snp_vcf[:-3] will be removed by bgzip.
    subprocess.run(["bgzip", "-f", out_snp_vcf[:-3]], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", out_snp_vcf], check=True)
logging.info("finished annotate_snps_pseudobulk")
