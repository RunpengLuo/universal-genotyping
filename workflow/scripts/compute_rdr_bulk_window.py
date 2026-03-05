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
from scipy.stats import pearsonr, gaussian_kde

from utils import *
from io_utils import read_region_file
from postprocess_utils import plot_1d_sample, get_mask_by_region_intervals

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
bias_bed_file = sm.input["bias_bed"]
genome_size = sm.input["genome_size"]
region_bed_file = sm.input["region_bed"]

sample_name = sm.params["sample_name"]
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
mosdepth_dir = sm.params["mosdepth_dir"]
chromosomes = sm.params["chromosomes"]

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
    repliseq=None,
    samplesize=50000,
    routlier=0.01,
    doutlier=0.001,
    min_mappability=0.9,
):
    """HMMcopy-style correctReadcount: GC + optional mappability + optional
    replication timing stages (ACEseq-style sequential LOWESS correction).

    Parameters
    ----------
    reads : np.ndarray
        Raw read counts per window (1-D).
    gc : np.ndarray
        Per-window GC fraction in [0, 1].
    mappability : np.ndarray or None
        Per-window mappability values in [0, 1]. If None, stage 2 is skipped.
    repliseq : np.ndarray or None
        Per-window consensus replication timing score. If None, stage 3 is
        skipped.  Higher values = earlier replication = higher expected
        coverage in cycling cells.
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
        ideal &= mappability >= min_mappability

    ideal_idx = np.where(ideal)[0]
    logging.info(
        f"  GC stage: {ideal_idx.size}/{int(valid.sum())} ideal bins ({ideal_idx.size / valid.sum() * 100:.1f}%)"
    )

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
        stage1[:, 0],
        stage1[:, 1],
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )
    grid_pred = stage1_interp(grid)
    stage2_smooth = lowess(grid_pred, grid, frac=0.3, return_sorted=True)
    final_gc_interp = interp1d(
        stage2_smooth[:, 0],
        stage2_smooth[:, 1],
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
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
    if mappability is not None:
        map_lo = np.nanquantile(mappability, doutlier)
        map_hi = np.nanquantile(mappability, 1.0 - doutlier)

        valid2 = (cor_gc > 0) & np.isfinite(mappability)
        read_hi2 = np.nanquantile(cor_gc[valid2], 1.0 - routlier)
        ideal2 = (
            valid2
            & (cor_gc <= read_hi2)
            & (mappability >= map_lo)
            & (mappability <= map_hi)
        )

        ideal_idx2 = np.where(ideal2)[0]
        logging.info(
            f"  MAP stage: {ideal_idx2.size}/{int(valid2.sum())} ideal bins ({ideal_idx2.size / valid2.sum() * 100:.1f}%)"
        )

        if ideal_idx2.size > samplesize:
            rng = np.random.default_rng(43)
            ideal_idx2 = rng.choice(ideal_idx2, size=samplesize, replace=False)

        map_ideal = mappability[ideal_idx2]
        cor_gc_ideal = cor_gc[ideal_idx2]

        s1_map = lowess(cor_gc_ideal, map_ideal, frac=0.03, return_sorted=True)
        grid_map = np.linspace(0, 1, 1001)
        s1_map_interp = interp1d(
            s1_map[:, 0],
            s1_map[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        grid_map_pred = s1_map_interp(grid_map)
        s2_map = lowess(grid_map_pred, grid_map, frac=0.3, return_sorted=True)
        final_map_interp = interp1d(
            s2_map[:, 0],
            s2_map[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
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

        cor_prev = cor_map
    else:
        cor_prev = cor_gc

    # --- Stage 3: Replication timing correction (ACEseq-style) ---
    if repliseq is None:
        return cor_prev.astype(np.float32)

    repli_finite = np.isfinite(repliseq)
    repli_lo = np.nanquantile(repliseq[repli_finite], doutlier)
    repli_hi = np.nanquantile(repliseq[repli_finite], 1.0 - doutlier)

    valid3 = (cor_prev > 0) & repli_finite
    read_hi3 = np.nanquantile(cor_prev[valid3], 1.0 - routlier)
    ideal3 = (
        valid3
        & (cor_prev <= read_hi3)
        & (repliseq >= repli_lo)
        & (repliseq <= repli_hi)
    )

    ideal_idx3 = np.where(ideal3)[0]
    logging.info(
        f"  REPLI stage: {ideal_idx3.size}/{int(valid3.sum())} ideal bins ({ideal_idx3.size / valid3.sum() * 100:.1f}%)"
    )

    if ideal_idx3.size > samplesize:
        rng = np.random.default_rng(44)
        ideal_idx3 = rng.choice(ideal_idx3, size=samplesize, replace=False)

    repli_ideal = repliseq[ideal_idx3]
    cor_prev_ideal = cor_prev[ideal_idx3]

    s1_repli = lowess(cor_prev_ideal, repli_ideal, frac=0.03, return_sorted=True)
    grid_repli = np.linspace(repli_lo, repli_hi, 1001)
    s1_repli_interp = interp1d(
        s1_repli[:, 0],
        s1_repli[:, 1],
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )
    grid_repli_pred = s1_repli_interp(grid_repli)
    s2_repli = lowess(grid_repli_pred, grid_repli, frac=0.3, return_sorted=True)
    final_repli_interp = interp1d(
        s2_repli[:, 0],
        s2_repli[:, 1],
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )

    repli_predicted = final_repli_interp(repliseq)
    with np.errstate(invalid="ignore", divide="ignore"):
        cor_repli = np.where(
            (repli_predicted > 0) & (cor_prev > 0),
            cor_prev / repli_predicted,
            0.0,
        )
    valid_repli = (repli_predicted > 0) & (cor_prev > 0)
    if valid_repli.any():
        scale3 = np.median(cor_prev[valid_repli]) / np.median(cor_repli[valid_repli])
        cor_repli[valid_repli] *= scale3

    return cor_repli.astype(np.float32)


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

    for title, dp_mat in [
        ("Before GC Correction", dp_before),
        ("After GC Correction", dp_after),
    ]:
        fig, axes = plt.subplots(
            1,
            nsamples,
            figsize=(panel_w, 5),
            squeeze=False,
        )
        for si, rep_id in enumerate(rep_ids):
            ax = axes[0, si]
            reads = dp_mat[:, si]
            valid = (reads > 0) & np.isfinite(gc)
            x, y = gc[valid], reads[valid]
            ylim = np.nanquantile(y, 0.99)

            # subsample for KDE fitting (speed), then evaluate on grid
            n_pts = len(x)
            rng = np.random.default_rng(0)
            if n_pts > 20000:
                idx = rng.choice(n_pts, size=20000, replace=False)
                kde = gaussian_kde(np.vstack([x[idx], y[idx]]))
            else:
                kde = gaussian_kde(np.vstack([x, y]))

            xgrid = np.linspace(x.min(), x.max(), 200)
            ygrid = np.linspace(0, ylim * 1.1, 200)
            xx, yy = np.meshgrid(xgrid, ygrid)
            positions = np.vstack([xx.ravel(), yy.ravel()])
            zz = kde(positions).reshape(xx.shape)

            ax.pcolormesh(xx, yy, zz, shading="gouraud", cmap="Blues", rasterized=True)
            ax.contour(xx, yy, zz, levels=6, colors="steelblue", linewidths=0.5, alpha=0.5)

            mad = np.median(np.abs(y - np.median(y)))
            r, _ = pearsonr(x, y)
            logging.info(f"  {title} {rep_id}: MAD={mad:.4f}  r={r:.4f}")
            ax.set_xlabel("GC Content")
            if si == 0:
                ax.set_ylabel("Observed Readcov")
            ax.set_title(f"{rep_id}\nMAD={mad:.4f}  r={r:.4f}")
            ax.set_xlim(x.min(), x.max())
            ax.set_ylim(0, ylim * 1.1)
        fig.suptitle(title, fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)


##################################################
logging.info("load GC BED")
bias_bed_full = pd.read_table(bias_bed_file, sep="\t")
assert "#CHR" in bias_bed_full.columns and "GC" in bias_bed_full.columns, (
    f"bias_bed must have #CHR, START, END, GC columns; got {bias_bed_full.columns.tolist()}"
)

##################################################
logging.info("load mosdepth depth data")
sample_df = pd.read_table(sample_file, sep="\t")
rep_ids = sample_df["REP_ID"].astype(str).tolist()
sample_types = sample_df["sample_type"].tolist()
nsamples = len(sample_df)
target_chroms = {f"chr{c}" for c in chromosomes}
join_keys = ["#CHR", "START", "END"]

mos_dfs = []
for rep_id in rep_ids:
    mos_file = os.path.join(mosdepth_dir, f"{rep_id}.regions.bed.gz")
    mos_df = pd.read_table(
        mos_file, sep="\t", header=None, names=["#CHR", "START", "END", "DEPTH"]
    )
    mos_df = mos_df[mos_df["#CHR"].isin(target_chroms)].reset_index(drop=True)
    mos_dfs.append(mos_df)

coords = mos_dfs[0][join_keys].copy()
n_windows = len(coords)
logging.info(f"{n_windows} windows across {len(target_chroms)} chromosomes")

bias_bed = pd.merge(left=coords, right=bias_bed_full, on=join_keys, how="left", sort=False)
n_matched = int(bias_bed["GC"].notna().sum())
logging.info(f"GC BED matched {n_matched}/{n_windows} windows")

# exclude centromeric windows (not in region_bed)
logging.info("building region mask from region_bed")
regions = read_region_file(region_bed_file)
win_mask_region = get_mask_by_region_intervals(bias_bed, regions)
n_excluded = int((~win_mask_region).sum())
logging.info(f"excluding {n_excluded}/{n_windows} centromeric windows")

bias_bed = bias_bed.loc[win_mask_region].reset_index(drop=True)
mos_dfs = [df.loc[win_mask_region].reset_index(drop=True) for df in mos_dfs]
n_windows = len(bias_bed)
logging.info(f"after filtering: {n_windows} windows")

gc_vals = bias_bed["GC"].to_numpy()
map_vals = bias_bed["MAP"].to_numpy() if "MAP" in bias_bed.columns else None
repli_vals = bias_bed["REPLI"].to_numpy(dtype=np.float64) if "REPLI" in bias_bed.columns else None
if repli_vals is not None:
    n_repli = int(np.isfinite(repli_vals).sum())
    logging.info(f"REPLI column found: {n_repli}/{n_windows} finite values")
else:
    logging.info("no REPLI column in bias_bed; skipping replication timing correction")

dp_mat = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, mos_df in enumerate(mos_dfs):
    dp_mat[:, i] = mos_df["DEPTH"].to_numpy(dtype=np.float32)
logging.info(f"depth matrix shape: {dp_mat.shape}")

##################################################
logging.info("applying correct_readcount per sample")
dp_corrected = np.zeros_like(dp_mat, dtype=np.float32)
for i, rep_id in enumerate(rep_ids):
    logging.info(f"correcting {rep_id}")
    dp_corrected[:, i] = correct_readcount(
        dp_mat[:, i],
        gc_vals,
        mappability=map_vals,
        repliseq=repli_vals,
        samplesize=samplesize,
        routlier=routlier,
        doutlier=doutlier,
        min_mappability=min_mappability,
    )

pdf = PdfPages(os.path.join(qc_dir, "window_depth_correction.pdf"))
plot_gc_correction_pdf(gc_vals, dp_mat, dp_corrected, rep_ids, pdf)
pdf.close()

rd_ylim = max(
    np.nanquantile(dp_mat, 0.99),
    np.nanquantile(dp_corrected, 0.99),
)
logging.info(f"depth plot shared y-lim={rd_ylim:.2f}")
for i, rep_id in enumerate(rep_ids):
    v = dp_mat[:, i]
    m = np.isfinite(v) & (v > 0)
    mad = np.median(np.abs(v[m] - np.median(v[m])))
    logging.info(f"  depth before correction {rep_id}: MAD={mad:.4f}")

    plot_file = os.path.join(qc_dir, f"depth_window_before_correction.{rep_id}.pdf")
    plot_1d_sample(
        bias_bed, v, genome_size, plot_file,
        unit="window", val_type="RD", max_ylim=rd_ylim,
    )

    vc = dp_corrected[:, i]
    mc = np.isfinite(vc) & (vc > 0)
    mad_c = np.median(np.abs(vc[mc] - np.median(vc[mc])))
    logging.info(f"  depth after correction {rep_id}: MAD={mad_c:.4f}")

    plot_file = os.path.join(qc_dir, f"depth_window_after_correction.{rep_id}.pdf")
    plot_1d_sample(
        bias_bed, vc, genome_size, plot_file,
        unit="window", val_type="RD", max_ylim=rd_ylim,
    )

np.savez_compressed(sm.output["dp_mtx"], mat=dp_corrected)

##################################################
has_normal = "normal" in sample_types
logging.info(f"compute RDRs, has_normal={has_normal}")

if has_normal:
    tumor_sidx = 1
    window_sizes = (bias_bed["END"] - bias_bed["START"]).to_numpy(dtype=np.float64)

    bases_mat = dp_corrected * window_sizes[:, None]
    total_bases = np.sum(bases_mat, axis=0)
    library_correction = total_bases[0] / total_bases[1:]
    logging.info(f"RDR library normalization factor: {library_correction}")

    normal_dp = dp_corrected[:, 0].copy()
    valid = np.isfinite(normal_dp) & (normal_dp > 0)
    logging.info(f"normal: {int(valid.sum())}/{n_windows} valid bins")

    with np.errstate(invalid="ignore", divide="ignore"):
        rdr_mat = dp_corrected[:, 1:] / normal_dp[:, None]
        rdr_mat *= library_correction[None, :]
    rdr_mat[~valid, :] = np.nan
else:
    tumor_sidx = 0
    # No matched normal: median-center each sample independently
    rdr_mat = np.full_like(dp_corrected, np.nan)
    for i in range(nsamples):
        col = dp_corrected[:, i]
        valid_i = np.isfinite(col) & (col > 0)
        med = np.median(col[valid_i])
        logging.info(f"  median-centering {rep_ids[i]}: median={med:.4f}")
        with np.errstate(invalid="ignore", divide="ignore"):
            rdr_mat[valid_i, i] = col[valid_i] / med

rdr_ylim = np.round(np.nanquantile(rdr_mat, 0.99)).astype(int) + 1
logging.info(f"window RDR y-lim={rdr_ylim}")

for i, rep_id in enumerate(rep_ids[tumor_sidx:]):
    v = rdr_mat[:, i]
    m = np.isfinite(v)
    mad = np.median(np.abs(v[m] - np.median(v[m])))
    logging.info(f"  window RDR {rep_id}: MAD={mad:.4f}")

    plot_file = os.path.join(qc_dir, f"rdr_window.{rep_id}.pdf")
    plot_1d_sample(
        bias_bed,
        v,
        genome_size,
        plot_file,
        unit="window",
        val_type="RDR",
        max_ylim=rdr_ylim,
    )

np.savez_compressed(sm.output["rdr_mtx"], mat=rdr_mat)

out_cols = ["#CHR", "START", "END", "GC"]
if "MAP" in bias_bed.columns:
    out_cols.append("MAP")
if "REPLI" in bias_bed.columns:
    out_cols.append("REPLI")
bias_bed[out_cols].to_csv(
    sm.output["bins_tsv"],
    sep="\t",
    header=True,
    index=False,
    compression="gzip",
)

##################################################
logging.info("aggregating window data into adaptive bins")
bbs = pd.read_table(bb_file, sep="\t")
n_bb = len(bbs)
n_tumors = rdr_mat.shape[1]
bb_rdr = np.full((n_bb, n_tumors), 1.0, dtype=np.float32)
bb_dp = np.full((n_bb, nsamples), np.nan, dtype=np.float32)

win_chr = bias_bed["#CHR"].to_numpy()
win_start = bias_bed["START"].to_numpy(dtype=np.int64)
win_end = bias_bed["END"].to_numpy(dtype=np.int64)

for chrom, bb_grp in bbs.groupby("#CHR", sort=False):
    win_mask = win_chr == chrom
    win_idx_chr = np.where(win_mask)[0]
    if len(win_idx_chr) == 0:
        continue

    bb_starts = bb_grp["START"].to_numpy(dtype=np.int64)
    bb_ends = bb_grp["END"].to_numpy(dtype=np.int64)
    bb_global_idx = bb_grp.index.to_numpy()
    n_bb_chr = len(bb_grp)

    win_s = win_start[win_idx_chr]
    win_e = win_end[win_idx_chr]

    # For each window, find all overlapping bins weighted by overlap fraction.
    # first_bin: leftmost bin overlapping the window (by window start)
    # last_bin:  rightmost bin overlapping the window (by window end - 1)
    first_bin = np.searchsorted(bb_ends, win_s, side="right")      # first bin whose end > win_start
    last_bin = np.searchsorted(bb_starts, win_e, side="left") - 1   # last bin whose start < win_end

    # Build (window_local_idx, bin_idx, overlap_weight) arrays
    ov_win_list, ov_bin_list, ov_wt_list = [], [], []
    for wi in range(len(win_idx_chr)):
        fb, lb = int(first_bin[wi]), int(last_bin[wi])
        if fb > lb or fb >= n_bb_chr or lb < 0:
            continue
        fb = max(fb, 0)
        lb = min(lb, n_bb_chr - 1)
        ws, we = int(win_s[wi]), int(win_e[wi])
        win_len = we - ws
        if win_len <= 0:
            continue
        for bi in range(fb, lb + 1):
            ov_start = max(ws, int(bb_starts[bi]))
            ov_end = min(we, int(bb_ends[bi]))
            ov_len = ov_end - ov_start
            if ov_len > 0:
                ov_win_list.append(wi)
                ov_bin_list.append(bi)
                ov_wt_list.append(ov_len / win_len)

    if len(ov_win_list) == 0:
        logging.info(f"  {chrom}: {n_bb_chr} bins, 0 windows assigned")
        continue

    ov_win = np.array(ov_win_list, dtype=np.intp)
    ov_bin = np.array(ov_bin_list, dtype=np.intp)
    ov_wt = np.array(ov_wt_list, dtype=np.float64)
    n_windows_assigned = len(np.unique(ov_win))

    # Aggregate corrected depth (all samples)
    for s in range(nsamples):
        dp_chr = dp_corrected[win_idx_chr, s]
        vals = dp_chr[ov_win]
        finite = np.isfinite(vals)
        bi_f = ov_bin[finite]
        wt_f = ov_wt[finite]
        v_f = vals[finite]
        sums = np.bincount(bi_f, weights=v_f * wt_f, minlength=n_bb_chr)
        wt_sums = np.bincount(bi_f, weights=wt_f, minlength=n_bb_chr)
        means = np.where(wt_sums > 0, sums / wt_sums, np.nan)
        bb_dp[bb_global_idx, s] = means

    # Aggregate RDR (tumor samples only)
    n_defaulted = 0
    for t in range(n_tumors):
        rdr_chr = rdr_mat[win_idx_chr, t]
        vals = rdr_chr[ov_win]
        finite = np.isfinite(vals)
        bi_f = ov_bin[finite]
        wt_f = ov_wt[finite]
        v_f = vals[finite]
        sums = np.bincount(bi_f, weights=v_f * wt_f, minlength=n_bb_chr)
        wt_sums = np.bincount(bi_f, weights=wt_f, minlength=n_bb_chr)
        means = np.where(wt_sums > 0, sums / wt_sums, 1.0)
        bb_rdr[bb_global_idx, t] = means
        if t == 0:
            n_defaulted = int((wt_sums == 0).sum())
    logging.info(
        f"  {chrom}: {n_bb_chr} bins, {n_windows_assigned} windows assigned, "
        f"{n_defaulted} bins defaulted to 1.0"
    )

n_filled = int(np.isfinite(bb_rdr[:, 0]).sum())
logging.info(f"bb RDR: {n_filled}/{n_bb} bins filled from window RDR")
np.savez_compressed(sm.output["rdr_mtx_bb"], mat=bb_rdr)
np.savez_compressed(sm.output["dp_mtx_bb"], mat=bb_dp)

for i, rep_id in enumerate(rep_ids[tumor_sidx:]):
    v = bb_rdr[:, i]
    m = np.isfinite(v)
    mad = np.median(np.abs(v[m] - np.median(v[m])))
    logging.info(f"  bb RDR {rep_id}: MAD={mad:.4f}")

    plot_file = os.path.join(qc_dir, f"rdr_bb.{rep_id}.pdf")
    plot_1d_sample(
        bbs,
        v,
        genome_size,
        plot_file,
        unit="bb",
        val_type="RDR",
        max_ylim=rdr_ylim,
    )

logging.info("finished.")
