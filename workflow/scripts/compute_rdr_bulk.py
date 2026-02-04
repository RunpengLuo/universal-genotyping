import os, sys, gzip, logging
from snakemake.script import snakemake as sm

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from utils import *
from bias_correction import *
from postprocess_utils import plot_1d_sample

def compute_RDR(
    bbs: pd.DataFrame,
    dp_mtx: np.ndarray,
    gc_dir: str,
    gc_correct=True,
    ref_file=None,
    mapp_file=None,
    genome_size=None,
    has_normal=True,
    tumor_sidx=1,
):
    assert has_normal, "no normal sample, TODO"
    bases_mat = dp_mtx * bbs["BLOCKSIZE"].to_numpy()[:, None]
    total_bases = np.sum(bases_mat, axis=0)
    library_correction = total_bases[0] / total_bases[1:]
    logging.info(f"RDR library normalization factor: {library_correction}")

    rdr_mat = dp_mtx[:, 1:] / dp_mtx[:, 0][:, None]
    rdr_mat *= library_correction[None, :]
    if gc_correct:
        gc_df = compute_gc_content(bbs, ref_file, mapp_file, genome_size)
        rdr_mat = bias_correction_rdr(rdr_mat, gc_df, gc_dir)
    return rdr_mat


##################################################
logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# input
sample_file = sm.input["sample_file"]
bb_file = sm.input["bb_file"]
reference = sm.input["reference"]
genome_size = sm.input["genome_size"]
mappability_file = maybe_path(sm.input["mappability_file"])

sample_name = sm.params["sample_name"]
qc_dir = sm.params["qc_dir"]
mosdepth_dir = sm.params["mosdepth_dir"]
gc_correct = bool(sm.params["gc_correct"])

logging.info("run compute_rdr_bulk")

##################################################
logging.info("concat mosdepth per-sample depth data")
samples_df = pd.read_table(sample_file, sep="\t")
rep_ids = samples_df["REP_ID"].astype(str).tolist()
nsamples = len(samples_df)
bbs = pd.read_table(bb_file, sep="\t")
for rep_id in rep_ids:
    mos_file = os.path.join(mosdepth_dir, f"{rep_id}.regions.bed.gz")
    mos_df = pd.read_table(
        mos_file, sep="\t", header=None, names=["#CHR", "START", "END", rep_id]
    )
    dup = mos_df.duplicated(subset=["#CHR", "START", "END"], keep=False)
    assert not dup.any(), f"{mos_file} has duplicate intervals: {int(dup.sum())}"
    assert len(mos_df) == len(bbs), (
        f"{mos_file} rowcount={len(mos_df)} != bbs rowcount={len(bbs)}"
    )
    bbs = bbs.merge(mos_df, how="left", on=["#CHR", "START", "END"], sort=False)
dp_mtx_bb = bbs[rep_ids].to_numpy(dtype=np.float32)
np.savez_compressed(sm.output["dp_mtx_bb"], mat=dp_mtx_bb)

##################################################
has_normal = "normal" in rep_ids
logging.info(f"compute RDRs, has_normal={has_normal}, gc_correct={gc_correct}")
assert has_normal, "no normal sample, TODO"
tumor_sidx = {False: 0, True: 1}[has_normal]

rdr_mtx_bb = compute_RDR(
    bbs,
    dp_mtx_bb,
    qc_dir,
    gc_correct=gc_correct,
    ref_file=reference,
    mapp_file=mappability_file,
    genome_size=genome_size,
    has_normal=has_normal,
    tumor_sidx=tumor_sidx,
)
np.savez_compressed(sm.output["rdr_mtx_bb"], mat=rdr_mtx_bb)

# plot per-sample RDRs
for i, rep_id in enumerate(rep_ids[tumor_sidx:]):
    plot_file = os.path.join(
        qc_dir, f"rdr_bb.{rep_id}.pdf"
    )
    plot_1d_sample(bbs, rdr_mtx_bb[:, i], genome_size, plot_file, unit="bb", val_type="RDR")
logging.info(f"finished.")
