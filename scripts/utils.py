import os
import sys
import subprocess
from io import StringIO

import pandas as pd
import numpy as np

from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy import sparse

def read_barcodes(bc_file: str):
    barcodes = []
    with open(bc_file, "r") as fd:
        for line in fd:
            barcodes.append(line.strip().split("\t")[0])
        fd.close()
    return barcodes

def read_VCF_cellsnp_err_header(vcf_file: str):
    """cellsnp-lite has issue with its header"""
    fields = "%CHROM\t%POS\t%INFO"
    names = ["#CHR", "POS", "INFO"]
    cmd = ["bcftools", "query", "-f", fields, vcf_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    snps = pd.read_csv(StringIO(result.stdout), sep="\t", header=None, names=names)

    def extract_info_field(info_str, key):
        for field in info_str.split(";"):
            if field.startswith(f"{key}="):
                return int(field.split("=")[1])
        return pd.NA

    snps["DP"] = snps["INFO"].apply(lambda x: extract_info_field(x, "DP"))
    snps = snps.drop(columns="INFO")
    return snps

# Load Allele count data
def load_cellsnp_files(
    cellsnp_dir: str,
    barcodes: list,
):
    print(f"load cell-snp files from {cellsnp_dir}")
    barcode_file = os.path.join(cellsnp_dir, "cellSNP.samples.tsv")
    vcf_file = os.path.join(cellsnp_dir, "cellSNP.base.vcf.gz")
    dp_file = os.path.join(cellsnp_dir, "cellSNP.tag.DP.mtx")
    ad_file = os.path.join(cellsnp_dir, "cellSNP.tag.AD.mtx")

    raw_barcodes = read_barcodes(barcode_file)  # assume no header
    barcode_indices = np.array([raw_barcodes.index(x) for x in barcodes])
    dp_mat: csr_matrix = mmread(dp_file).tocsr()
    alt_mat: csr_matrix = mmread(ad_file).tocsr()
    ref_mat = dp_mat - alt_mat

    dp_mat = dp_mat[:, barcode_indices]
    alt_mat = alt_mat[:, barcode_indices]
    ref_mat = ref_mat[:, barcode_indices]

    cell_snps = read_VCF_cellsnp_err_header(vcf_file)
    cell_snps["RAW_SNP_IDX"] = np.arange(len(cell_snps))  # use to index matrix
    return [cell_snps, dp_mat, ref_mat, alt_mat]
