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
from scipy.interpolate import interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import pearsonr

from utils import *
from postprocess_utils import plot_1d_sample

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

##################################################
logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# input
sample_file = sm.input["sample_file"]
bb_file = sm.input["bb_file"]
gc_bed_file = sm.input["gc_bed"]
genome_size = sm.input["genome_size"]
mappability_file = maybe_path(sm.input["mappability_file"])

sample_name = sm.params["sample_name"]
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
mosdepth_dir = sm.params["mosdepth_dir"]

# params
samplesize = int(sm.params["samplesize"])
routlier = float(sm.params["routlier"])
doutlier = float(sm.params["doutlier"])
min_mappability = float(sm.params["min_mappability"])

logging.info("run compute_rdr_bulk_window")


##################################################
def correct_readcount(
    reads,
    gc,
    mappability=None,
    samplesize=50000,
    routlier=0.01,
    doutlier=0.001,
    min_mappability=0.9,
):
    """Full HMMcopy-style correctReadcount: GC stage + optional mappability stage.

    Parameters
    ----------
    reads : np.ndarray
        Raw read counts per window (1-D).
    gc : np.ndarray
        Per-window GC fraction in [0, 1].
    mappability : np.ndarray or None
        Per-window mappability values in [0, 1]. If None, stage 2 is skipped.
    samplesize : int
        Max number of ideal bins for loess fitting.
    routlier : float
        Upper quantile for read-count outlier removal.
    doutlier : float
        Quantile for GC/mappability domain outlier removal (top/bottom).
    min_mappability : float
        Minimum mappability for ideal bins.

    Returns
    -------
    np.ndarray
        Corrected read counts (same length as *reads*).
    """
    reads = reads.astype(np.float64)
    n = len(reads)

    # --- Stage 1: GC correction ---
    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)

    valid = (reads > 0) & np.isfinite(gc) & (gc >= 0)
    read_hi = np.nanquantile(reads[valid], 1.0 - routlier)
    ideal = valid & (reads <= read_hi) & (gc >= gc_lo) & (gc <= gc_hi)
    if mappability is not None:
        ideal &= (mappability >= min_mappability)

    ideal_idx = np.where(ideal)[0]
    logging.info(f"  GC stage: {ideal_idx.size} ideal bins out of {int(valid.sum())} valid")

    if ideal_idx.size > samplesize:
        rng = np.random.default_rng(42)
        ideal_idx = rng.choice(ideal_idx, size=samplesize, replace=False)

    gc_ideal = gc[ideal_idx]
    reads_ideal = reads[ideal_idx]

    # Stage 1a: tight loess on ideal bins
    stage1 = lowess(reads_ideal, gc_ideal, frac=0.03, return_sorted=True)
    # Stage 1b: smooth on fine grid
    grid = np.linspace(0, 1, 1001)
    stage1_interp = interp1d(
        stage1[:, 0], stage1[:, 1],
        kind="linear", bounds_error=False, fill_value="extrapolate",
    )
    grid_pred = stage1_interp(grid)
    stage2_smooth = lowess(grid_pred, grid, frac=0.3, return_sorted=True)
    final_gc_interp = interp1d(
        stage2_smooth[:, 0], stage2_smooth[:, 1],
        kind="linear", bounds_error=False, fill_value="extrapolate",
    )

    gc_predicted = final_gc_interp(gc)
    with np.errstate(invalid="ignore", divide="ignore"):
        cor_gc = np.where(
            (gc_predicted > 0) & (reads > 0),
            reads / gc_predicted,
            0.0,
        )
    # Rescale so corrected median matches original median
    valid_gc = (gc_predicted > 0) & (reads > 0)
    if valid_gc.any():
        scale = np.median(reads[valid_gc]) / np.median(cor_gc[valid_gc])
        cor_gc[valid_gc] *= scale

    # --- Stage 2: Mappability correction ---
    if mappability is None:
        return cor_gc.astype(np.float32)

    map_lo = np.nanquantile(mappability, doutlier)
    map_hi = np.nanquantile(mappability, 1.0 - doutlier)

    valid2 = (cor_gc > 0) & np.isfinite(mappability)
    read_hi2 = np.nanquantile(cor_gc[valid2], 1.0 - routlier)
    ideal2 = valid2 & (cor_gc <= read_hi2) & (mappability >= map_lo) & (mappability <= map_hi)

    ideal_idx2 = np.where(ideal2)[0]
    logging.info(f"  MAP stage: {ideal_idx2.size} ideal bins out of {int(valid2.sum())} valid")

    if ideal_idx2.size > samplesize:
        rng = np.random.default_rng(43)
        ideal_idx2 = rng.choice(ideal_idx2, size=samplesize, replace=False)

    map_ideal = mappability[ideal_idx2]
    cor_gc_ideal = cor_gc[ideal_idx2]

    s1_map = lowess(cor_gc_ideal, map_ideal, frac=0.03, return_sorted=True)
    grid_map = np.linspace(0, 1, 1001)
    s1_map_interp = interp1d(
        s1_map[:, 0], s1_map[:, 1],
        kind="linear", bounds_error=False, fill_value="extrapolate",
    )
    grid_map_pred = s1_map_interp(grid_map)
    s2_map = lowess(grid_map_pred, grid_map, frac=0.3, return_sorted=True)
    final_map_interp = interp1d(
        s2_map[:, 0], s2_map[:, 1],
        kind="linear", bounds_error=False, fill_value="extrapolate",
    )

    map_predicted = final_map_interp(mappability)
    with np.errstate(invalid="ignore", divide="ignore"):
        cor_map = np.where(
            (map_predicted > 0) & (cor_gc > 0),
            cor_gc / map_predicted,
            0.0,
        )
    valid_map = (map_predicted > 0) & (cor_gc > 0)
    if valid_map.any():
        scale2 = np.median(cor_gc[valid_map]) / np.median(cor_map[valid_map])
        cor_map[valid_map] *= scale2

    return cor_map.astype(np.float32)


def plot_gc_correction_pdf(gc, dp_before, dp_after, rep_ids, pdf):
    """Two-page PDF: page 1 = before GC correction, page 2 = after.

    Each page has one subplot per sample (side by side), showing a
    density-style scatter of GC Content vs Observed Readcov.

    Parameters
    ----------
    gc : np.ndarray
        Per-window GC fraction.
    dp_before : np.ndarray
        Raw depth matrix (windows x samples).
    dp_after : np.ndarray
        Corrected depth matrix (windows x samples).
    rep_ids : list[str]
        Sample/replicate identifiers.
    pdf : PdfPages
        Open PdfPages handle.
    """
    nsamples = len(rep_ids)
    panel_w = max(5, 5 * nsamples)

    for title, dp_mat in [("Before GC Correction", dp_before),
                          ("After GC Correction", dp_after)]:
        fig, axes = plt.subplots(
            1, nsamples, figsize=(panel_w, 5), squeeze=False,
        )
        for si, rep_id in enumerate(rep_ids):
            ax = axes[0, si]
            reads = dp_mat[:, si]
            valid = (reads > 0) & np.isfinite(gc)
            ax.scatter(
                gc[valid], reads[valid], s=1, alpha=0.05,
                rasterized=True, color="steelblue",
            )
            mad = np.median(np.abs(reads[valid] - np.median(reads[valid])))
            r, _ = pearsonr(gc[valid], reads[valid])
            ax.set_xlabel("GC Content")
            if si == 0:
                ax.set_ylabel("Observed Readcov")
            ax.set_title(f"{rep_id}\nMAD={mad:.4f}  r={r:.4f}")
            ax.set_xlim(0, 1)
        fig.suptitle(title, fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)


##################################################
logging.info("load GC BED")
gc_bed = pd.read_table(gc_bed_file, sep="\t")
assert "#CHR" in gc_bed.columns and "GC" in gc_bed.columns, (
    f"gc_bed must have #CHR, START, END, GC columns; got {gc_bed.columns.tolist()}"
)
gc_vals = gc_bed["GC"].to_numpy()
map_vals = gc_bed["MAP"].to_numpy() if "MAP" in gc_bed.columns else None

# override with external mappability if provided
if mappability_file is not None:
    logging.info(f"loading external mappability from {mappability_file}")
    ext_map = pd.read_table(mappability_file, sep="\t")
    if "MAP" in ext_map.columns:
        map_merged = gc_bed[["#CHR", "START", "END"]].merge(
            ext_map.rename(columns={ext_map.columns[3]: "MAP"})[
                [ext_map.columns[0], ext_map.columns[1], ext_map.columns[2], "MAP"]
            ].rename(columns={ext_map.columns[0]: "#CHR", ext_map.columns[1]: "START", ext_map.columns[2]: "END"}),
            on=["#CHR", "START", "END"], how="left",
        )
        map_vals = pd.to_numeric(map_merged["MAP"], errors="coerce").fillna(1.0).clip(0.0, 1.0).to_numpy()
    else:
        logging.warning("external mappability file has no MAP column; using gc_bed MAP")

n_windows = len(gc_bed)
logging.info(f"GC BED: {n_windows} windows")

##################################################
logging.info("concat mosdepth per-sample depth data")
sample_df = pd.read_table(sample_file, sep="\t")
rep_ids = sample_df["REP_ID"].astype(str).tolist()
sample_types = sample_df["sample_type"].tolist()
nsamples = len(sample_df)

dp_mat = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, rep_id in enumerate(rep_ids):
    mos_file = os.path.join(mosdepth_dir, f"{rep_id}.regions.bed.gz")
    mos_df = pd.read_table(
        mos_file, sep="\t", header=None, names=["#CHR", "START", "END", "DEPTH"]
    )
    assert len(mos_df) == n_windows, (
        f"{mos_file} rowcount={len(mos_df)} != gc_bed rowcount={n_windows}"
    )
    dp_mat[:, i] = mos_df["DEPTH"].to_numpy(dtype=np.float32)

logging.info(f"depth matrix shape: {dp_mat.shape}")

##################################################
# Apply HMMcopy correctReadcount per sample
logging.info("applying correct_readcount per sample")
dp_corrected = np.zeros_like(dp_mat, dtype=np.float32)
for i, rep_id in enumerate(rep_ids):
    logging.info(f"correcting {rep_id}")
    dp_corrected[:, i] = correct_readcount(
        dp_mat[:, i], gc_vals,
        mappability=map_vals,
        samplesize=samplesize,
        routlier=routlier,
        doutlier=doutlier,
        min_mappability=min_mappability,
    )

# QC: before/after GC correction (page 1 = before, page 2 = after)
pdf = PdfPages(os.path.join(qc_dir, "window_depth_correction.pdf"))
plot_gc_correction_pdf(gc_vals, dp_mat, dp_corrected, rep_ids, pdf)
pdf.close()

np.savez_compressed(sm.output["dp_mtx"], mat=dp_corrected)

##################################################
# Compute RDR: tumor/normal ratio with library-size normalization
has_normal = "normal" in sample_types
logging.info(f"compute RDRs, has_normal={has_normal}")
assert has_normal, "no normal sample for RDR computation"
tumor_sidx = 1  # normal is index 0

window_sizes = (gc_bed["END"] - gc_bed["START"]).to_numpy(dtype=np.float64)

bases_mat = dp_corrected * window_sizes[:, None]
total_bases = np.sum(bases_mat, axis=0)
library_correction = total_bases[0] / total_bases[1:]
logging.info(f"RDR library normalization factor: {library_correction}")

normal_dp = dp_corrected[:, 0].copy()
mask = np.isfinite(normal_dp) & (normal_dp > 0)
num_valid = int(mask.sum())
if num_valid < n_windows:
    logging.warning(f"#invalid normal bins={n_windows - num_valid}/{n_windows}")
    if num_valid == 0:
        raise ValueError("All normal bins invalid; cannot compute RDR.")
    fill = float(np.nanmedian(normal_dp[mask]))
    logging.warning(f"impute normal depth via nanmedian: {fill}")
    normal_dp[~mask] = fill

rdr_mat = dp_corrected[:, 1:] / normal_dp[:, None]
rdr_mat *= library_correction[None, :]

rdr_ylim = np.round(np.nanquantile(rdr_mat, 0.99)).astype(int) + 1
logging.info(f"RDR y-lim={rdr_ylim}")

# QC plots: per-sample RDR
for i, rep_id in enumerate(rep_ids[tumor_sidx:]):
    plot_file = os.path.join(qc_dir, f"rdr_window.{rep_id}.pdf")
    plot_1d_sample(
        gc_bed,
        rdr_mat[:, i],
        genome_size,
        plot_file,
        unit="window",
        val_type="RDR",
        max_ylim=rdr_ylim,
    )

np.savez_compressed(sm.output["rdr_mtx"], mat=rdr_mat)

# Write bin info
gc_bed[["#CHR", "START", "END", "GC", "MAP"]].to_csv(
    sm.output["bins_tsv"], sep="\t", header=True, index=False,
    compression="gzip",
)

##################################################
# Aggregate window-level RDR into adaptive allele bins (bb)
logging.info("aggregating window RDR into adaptive bins")
bbs = pd.read_table(bb_file, sep="\t")
n_bb = len(bbs)
n_tumors = rdr_mat.shape[1]
bb_rdr = np.full((n_bb, n_tumors), np.nan, dtype=np.float32)

win_chr = gc_bed["#CHR"].to_numpy()
win_start = gc_bed["START"].to_numpy(dtype=np.int64)

for chrom, bb_grp in bbs.groupby("#CHR", sort=False):
    win_mask = win_chr == chrom
    win_idx_chr = np.where(win_mask)[0]
    if len(win_idx_chr) == 0:
        continue

    bb_starts = bb_grp["START"].to_numpy(dtype=np.int64)
    bb_ends = bb_grp["END"].to_numpy(dtype=np.int64)
    bb_global_idx = bb_grp.index.to_numpy()
    n_bb_chr = len(bb_grp)

    # Assign each window to a bin via vectorised searchsorted
    win_starts_chr = win_start[win_idx_chr]
    bin_idx = np.searchsorted(bb_starts, win_starts_chr, side="right") - 1
    valid_assign = (bin_idx >= 0) & (bin_idx < n_bb_chr)
    valid_assign[valid_assign] &= win_starts_chr[valid_assign] < bb_ends[bin_idx[valid_assign]]

    for t in range(n_tumors):
        rdr_chr = rdr_mat[win_idx_chr, t]
        finite = valid_assign & np.isfinite(rdr_chr)
        bi = bin_idx[finite]
        vals = rdr_chr[finite]
        sums = np.bincount(bi, weights=vals, minlength=n_bb_chr)
        counts = np.bincount(bi, minlength=n_bb_chr)
        means = np.where(counts > 0, sums / counts, np.nan)
        bb_rdr[bb_global_idx, t] = means

n_filled = int(np.isfinite(bb_rdr[:, 0]).sum())
logging.info(f"bb RDR: {n_filled}/{n_bb} bins filled from window RDR")
np.savez_compressed(sm.output["rdr_mtx_bb"], mat=bb_rdr)

logging.info("finished.")
