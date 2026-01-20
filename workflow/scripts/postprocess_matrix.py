import os, sys, gzip, logging
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

logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

logging.info("postprocess count matrix")

data_types = sm.params["data_types"]
rep_ids = sm.params["rep_ids"]

phased_snps = read_VCF(sm.input["snp_file"])
phased_snps["KEY"] = (
    phased_snps["CHROM"].astype(str) + "_" + phased_snps["POS"].astype(str)
)

gt = phased_snps["GT"].astype(str)
is_phased = gt.str.contains(r"\|", na=False).to_numpy()
phases = None
if is_phased.any():
    logging.info("SNP file is phased")
    phased_snps = phased_snps.loc[is_phased].reset_index(drop=True)
    phases = phased_snps["GT"].astype(str).str[2].astype(np.int8).to_numpy()

parent_keys = pd.Index(phased_snps["KEY"])
assert not parent_keys.duplicated().any()

M = len(parent_keys)
logging.info(f"total #SNPs={M}")

for data_type in data_types:
    barcode_list = []
    dp_list = []
    ad_list = []
    for rep_id in rep_ids:
        barcodes = pd.read_table(
            f"pileup/{data_type}_{rep_id}/cellSNP.samples.tsv",
            sep="\t",
            header=None,
            names=["BARCODE"],
        )
        if len(rep_ids) > 1:
            barcodes["BARCODE"] = barcodes["BARCODE"].astype(str) + f"_{rep_id}"
        child_snps = read_VCF(f"pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz")
        m = len(child_snps)
        logging.info(f"MOD={data_type}, REP_ID={rep_id}, #SNPs={m}")
        child_snps["KEY"] = (
            child_snps["CHROM"].astype(str) + "_" + child_snps["POS"].astype(str)
        )
        child_keys = pd.Index(child_snps["KEY"])
        assert not child_keys.duplicated().any()

        # child slice SNP location in parent, -1 if not presented in child (M,)
        # e.g., [0, *, 1, 2, *], M=5, m=3
        child_loc = child_keys.get_indexer(parent_keys)

        dp_mtx = scipy.io.mmread(f"pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx").tocsr()
        ad_mtx = scipy.io.mmread(f"pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx").tocsr()

        assert dp_mtx.shape == ad_mtx.shape
        assert dp_mtx.shape[0] == len(child_snps)
        assert dp_mtx.shape[1] == len(barcodes)

        present = child_loc >= 0
        logging.info(
            f"matched SNPs in parent={present.sum()}/{M} (child #SNPs={len(child_keys)})"
        )
        n_cells = dp_mtx.shape[1]

        # parent-row indices that exist in child slice
        row_new = np.flatnonzero(present)  # indices in [0..M),
        row_old = child_loc[present]  # corresponding indices in child [0..n_child)

        # DP
        dp_present = dp_mtx[row_old, :].tocoo()
        dp_canon = scipy.sparse.csr_matrix(
            (dp_present.data, (row_new[dp_present.row], dp_present.col)),
            shape=(M, n_cells),
            dtype=dp_mtx.dtype,
        )

        # AD
        ad_present = ad_mtx[row_old, :].tocoo()
        ad_canon = scipy.sparse.csr_matrix(
            (ad_present.data, (row_new[ad_present.row], ad_present.col)),
            shape=(M, n_cells),
            dtype=ad_mtx.dtype,
        )

        barcode_list.append(barcodes)
        dp_list.append(dp_canon)
        ad_list.append(ad_canon)

    # merge replicates
    all_barcodes = pd.concat(barcode_list, axis=0, ignore_index=True)
    dp_mtx = scipy.sparse.hstack(dp_list, format="csr")
    alt_mtx = scipy.sparse.hstack(ad_list, format="csr")
    ref_mtx = dp_mtx - alt_mtx
    assert (ref_mtx.data >= 0).all(), "invalid data, DP < AD?"
    if phases is None:
        # no phasing, default ALT=A, REF=B
        a_mtx = alt_mtx
        b_mtx = ref_mtx
    else:
        b_mtx = ref_mtx.multiply(phases[:, None]) + alt_mtx.multiply(
            1 - phases[:, None]
        )
        a_mtx = dp_mtx - b_mtx

    row_dp = np.asarray(dp_mtx.sum(axis=1)).ravel()
    a_mtx = a_mtx[row_dp > 0, :]
    b_mtx = b_mtx[row_dp > 0, :]
    snp_ids = parent_keys[row_dp > 0].to_numpy()

    dp_mtx = dp_mtx[row_dp > 0, :]
    num_snps = dp_mtx.shape[0]
    num_barcodes = len(all_barcodes)
    sparsity = 1.0 - dp_mtx.nnz / (num_barcodes * num_snps)
    logging.info(
        f"final matrix: #barcodes={num_barcodes}, #SNPs={num_snps}, sparsity={sparsity:.3f}"
    )

    all_barcodes.to_csv(f"{data_type}/barcodes.txt", header=False, index=False)
    scipy.sparse.save_npz(f"{data_type}/cell_snp_Aallele.npz", a_mtx)
    scipy.sparse.save_npz(f"{data_type}/cell_snp_Ballele.npz", b_mtx)
    np.save(f"{data_type}/unique_snp_ids.npy", snp_ids)

logging.info("finished postprocess_matrix")
