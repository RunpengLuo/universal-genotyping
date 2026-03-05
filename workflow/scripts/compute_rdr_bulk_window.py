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
from postprocess_utils import (
    plot_1d_sample,
    get_mask_by_region_intervals,
    compute_overlap_weights,
    weighted_bincount_mean,
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

logging.basicConfig(
    filename=sm.log[0],
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def correct_readcount(
    reads,
    gc,
    mappability=None,
    repliseq=None,
    samplesize=50000,
    routlier=0.01,
    doutlier=0.001,
    min_mappability=0.9,
    lowess_frac_tight=0.03,
    lowess_frac_smooth=0.3,
    grid_size=1001,
    seed=42,
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
    lowess_frac_tight : float
        Bandwidth fraction for the first (tight) LOWESS pass.
    lowess_frac_smooth : float
        Bandwidth fraction for the second (smoothing) LOWESS pass.
    grid_size : int
        Number of grid points for LOWESS interpolation.
    seed : int
        Base random seed for subsampling ideal bins (incremented per stage).

    Returns
    -------
    np.ndarray
        Corrected read counts (same length as *reads*).  NaN where
        correction is not possible.
    """
    reads = reads.astype(np.float64)
    n = len(reads)

    def _fit_lowess_interp(y, x, grid):
        """Tight LOWESS -> grid smooth -> final interpolator."""
        s1 = lowess(y, x, frac=lowess_frac_tight, return_sorted=True)
        s1_fn = interp1d(
            s1[:, 0],
            s1[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        s2 = lowess(s1_fn(grid), grid, frac=lowess_frac_smooth, return_sorted=True)
        return interp1d(
            s2[:, 0],
            s2[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )

    def _apply_stage(
        prev, covariate, cov_lo, cov_hi, grid, seed,
        stage_name="", extra_ideal_mask=None,
    ):
        """Apply one LOWESS correction stage and return corrected values."""
        valid = (prev > 0) & np.isfinite(covariate)
        val_hi = np.nanquantile(prev[valid], 1.0 - routlier)
        ideal = valid & (prev <= val_hi) & (covariate >= cov_lo) & (covariate <= cov_hi)
        if extra_ideal_mask is not None:
            ideal &= extra_ideal_mask

        ideal_idx = np.where(ideal)[0]
        if ideal_idx.size > samplesize:
            rng = np.random.default_rng(seed)
            ideal_idx = rng.choice(ideal_idx, size=samplesize, replace=False)

        interp_fn = _fit_lowess_interp(prev[ideal_idx], covariate[ideal_idx], grid)
        predicted = interp_fn(covariate)

        with np.errstate(invalid="ignore", divide="ignore"):
            corrected = np.where(
                (predicted > 0) & (prev > 0),
                prev / predicted,
                np.nan,
            )
        valid_corr = (predicted > 0) & (prev > 0)
        if valid_corr.any():
            scale = np.median(prev[valid_corr]) / np.median(corrected[valid_corr])
            corrected[valid_corr] *= scale

        n_valid = int(valid.sum())
        n_ideal = ideal_idx.size
        n_nan = int(np.isnan(corrected).sum())
        ideal_pct = n_ideal / max(n_valid, 1) * 100
        nan_pct = n_nan / max(n, 1) * 100
        logging.info(
            f"    {stage_name:<5s}  {n_ideal:>8d}/{n_valid} ({ideal_pct:5.1f}%) fit,  "
            f"{n_nan:>8d}/{n} ({nan_pct:5.1f}%) NaN"
        )
        return corrected

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)
    grid_01 = np.linspace(0, 1, grid_size)
    extra = (mappability >= min_mappability) if mappability is not None else None
    cor = _apply_stage(
        reads, gc, gc_lo, gc_hi, grid_01, seed=seed,
        stage_name="GC", extra_ideal_mask=extra,
    )

    if mappability is not None:
        map_lo = np.nanquantile(mappability, doutlier)
        map_hi = np.nanquantile(mappability, 1.0 - doutlier)
        cor = _apply_stage(
            cor, mappability, map_lo, map_hi, grid_01, seed=seed + 1,
            stage_name="MAP",
        )

    if repliseq is not None:
        repli_finite = np.isfinite(repliseq)
        repli_lo = np.nanquantile(repliseq[repli_finite], doutlier)
        repli_hi = np.nanquantile(repliseq[repli_finite], 1.0 - doutlier)
        grid_repli = np.linspace(repli_lo, repli_hi, grid_size)
        cor = _apply_stage(
            cor, repliseq, repli_lo, repli_hi, grid_repli, seed=seed + 2,
            stage_name="REPLI",
        )

    return cor.astype(np.float32)


def plot_gc_correction_pdf(gc, dp_before, dp_after, rep_ids, pdf):
    """Two-page PDF: before/after GC correction KDE density plots."""
    nsamples = len(rep_ids)
    panel_w = max(5, 5 * nsamples)

    for title, dp_mat in [
        ("Before GC Correction", dp_before),
        ("After GC Correction", dp_after),
    ]:
        fig, axes = plt.subplots(1, nsamples, figsize=(panel_w, 5), squeeze=False)
        for si, rep_id in enumerate(rep_ids):
            ax = axes[0, si]
            reads = dp_mat[:, si]
            valid = (reads > 0) & np.isfinite(gc)
            x, y = gc[valid], reads[valid]
            ylim = np.nanquantile(y, 0.99)

            n_pts = len(x)
            rng = np.random.default_rng(0)
            if n_pts > 20000:
                idx = rng.choice(n_pts, size=20000, replace=False)
                kde = gaussian_kde(np.vstack([x[idx], y[idx]]))
            else:
                kde = gaussian_kde(np.vstack([x, y]))

            xgrid = np.linspace(x.min(), x.max(), 200)
            ygrid = np.linspace(0, ylim * 1.5, 200)
            xx, yy = np.meshgrid(xgrid, ygrid)
            positions = np.vstack([xx.ravel(), yy.ravel()])
            zz = kde(positions).reshape(xx.shape)

            ax.pcolormesh(xx, yy, zz, shading="gouraud", cmap="Blues", rasterized=True)
            ax.contour(
                xx, yy, zz, levels=6, colors="steelblue", linewidths=0.5, alpha=0.5
            )

            mad = np.median(np.abs(y - np.median(y)))
            r, _ = pearsonr(x, y)
            logging.info(f"  {title:<24s} {rep_id:<8s}: MAD={mad:.4f}  r={r:.4f}")
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


def log_nan_summary(name, mat, labels, n_total):
    """Log per-column NaN counts for a matrix."""
    for i, label in enumerate(labels):
        col = mat[:, i] if mat.ndim == 2 else mat
        n_nan = int(np.isnan(col).sum())
        pct = n_nan / max(n_total, 1) * 100
        logging.info(f"  {name:<16s} {label:<8s}: {n_nan:>8d}/{n_total} ({pct:5.1f}%) NaN")


def log_mad_and_plot(
    pos_df,
    mat,
    labels,
    genome_size,
    qc_dir,
    prefix,
    unit,
    val_type,
    ylim,
    **plot_kwargs,
):
    """Log MAD per column and generate 1-D chromosome scatter plots."""
    for i, label in enumerate(labels):
        v = mat[:, i]
        m = np.isfinite(v)
        if m.any():
            mad = np.median(np.abs(v[m] - np.median(v[m])))
            logging.info(f"  {prefix:<28s} {label:<8s}: MAD={mad:.4f}")
        plot_file = os.path.join(qc_dir, f"{prefix}.{label}.pdf")
        plot_1d_sample(
            pos_df,
            v,
            genome_size,
            plot_file,
            unit=unit,
            val_type=val_type,
            max_ylim=ylim,
            **plot_kwargs,
        )


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

logging.info("run compute_rdr_bulk_window")

sample_file = sm.input["sample_file"]
bb_file = sm.input["bb_file"]
bias_bed_file = sm.input["bias_bed"]
genome_size = sm.input["genome_size"]
region_bed_file = sm.input["region_bed"]
blacklist_bed_file = maybe_path(sm.input["blacklist_bed"])

sample_name = sm.params["sample_name"]
qc_dir = sm.output["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
win_qc_dir = os.path.join(qc_dir, "window")
bb_qc_dir = os.path.join(qc_dir, "bb")
os.makedirs(win_qc_dir, exist_ok=True)
os.makedirs(bb_qc_dir, exist_ok=True)
mosdepth_dir = sm.params["mosdepth_dir"]
chromosomes = sm.params["chromosomes"]

samplesize = int(sm.params["samplesize"])
routlier = float(sm.params["routlier"])
doutlier = float(sm.params["doutlier"])
min_mappability = float(sm.params["min_mappability"])
gc_correct = bool(sm.params["gc_correct"])
rt_correct = bool(sm.params["rt_correct"])

# --- Load windows ---
logging.info("load bias BED and mosdepth depth")
bias_bed_full = pd.read_table(bias_bed_file, sep="\t")
assert "#CHR" in bias_bed_full.columns and "GC" in bias_bed_full.columns, (
    f"bias_bed must have #CHR, START, END, GC columns; got {bias_bed_full.columns.tolist()}"
)

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

bias_bed = pd.merge(
    left=coords, right=bias_bed_full, on=join_keys, how="left", sort=False
)
_gc_matched = int(bias_bed["GC"].notna().sum())
logging.info(f"GC BED matched {_gc_matched}/{n_windows} ({_gc_matched / max(n_windows, 1) * 100:.1f}%)")

dp_raw = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, mos_df in enumerate(mos_dfs):
    dp_raw[:, i] = mos_df["DEPTH"].to_numpy(dtype=np.float32)
rd_raw_ylim = max(np.nanquantile(dp_raw, 0.99), 1.0)
log_mad_and_plot(
    bias_bed,
    dp_raw,
    rep_ids,
    genome_size,
    win_qc_dir,
    "depth_before_correction",
    "window",
    "RD",
    rd_raw_ylim,
    region_bed=region_bed_file,
    blacklist_bed=blacklist_bed_file,
)

# --- Filter windows ---
regions = read_region_file(region_bed_file)
win_mask_region = get_mask_by_region_intervals(bias_bed, regions)
_n_excl = int((~win_mask_region).sum())
logging.info(
    f"region filter: excluding {_n_excl}/{n_windows} ({_n_excl / max(n_windows, 1) * 100:.1f}%)"
)

bias_bed = bias_bed.loc[win_mask_region].reset_index(drop=True)
mos_dfs = [df.loc[win_mask_region].reset_index(drop=True) for df in mos_dfs]
n_windows = len(bias_bed)

if blacklist_bed_file is not None:
    bl_regions = read_region_file(blacklist_bed_file)
    bl_mask = get_mask_by_region_intervals(bias_bed, bl_regions)
    n_bl = int(bl_mask.sum())
    logging.info(f"blacklist filter: excluding {n_bl}/{n_windows} ({n_bl / max(n_windows, 1) * 100:.1f}%)")
    bias_bed = bias_bed.loc[~bl_mask].reset_index(drop=True)
    mos_dfs = [df.loc[~bl_mask].reset_index(drop=True) for df in mos_dfs]
    n_windows = len(bias_bed)

logging.info(f"{n_windows} windows after filtering")

gc_vals = bias_bed["GC"].to_numpy()
map_vals = bias_bed["MAP"].to_numpy() if gc_correct and "MAP" in bias_bed.columns else None
repli_vals = (
    bias_bed["REPLI"].to_numpy(dtype=np.float64)
    if rt_correct and "REPLI" in bias_bed.columns
    else None
)
if repli_vals is not None:
    _n_finite = int(np.isfinite(repli_vals).sum())
    logging.info(
        f"REPLI column: {_n_finite}/{n_windows} ({_n_finite / max(n_windows, 1) * 100:.1f}%) finite"
    )
else:
    logging.info("no REPLI column; skipping replication timing correction")

dp_mat = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, mos_df in enumerate(mos_dfs):
    dp_mat[:, i] = mos_df["DEPTH"].to_numpy(dtype=np.float32)
logging.info(f"depth matrix shape: {dp_mat.shape}")

# --- Bias correction ---
if gc_correct:
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
else:
    logging.info("gc_correct=False; skipping bias correction")
    dp_corrected = dp_mat.copy()

log_nan_summary("corrected depth", dp_corrected, rep_ids, n_windows)

pdf = PdfPages(os.path.join(win_qc_dir, "depth_correction.pdf"))
plot_gc_correction_pdf(gc_vals, dp_mat, dp_corrected, rep_ids, pdf)
pdf.close()

rd_ylim = max(np.nanquantile(dp_corrected, 0.99), 1.0)
log_mad_and_plot(
    bias_bed,
    dp_corrected,
    rep_ids,
    genome_size,
    win_qc_dir,
    "depth_after_correction",
    "window",
    "RD",
    rd_ylim,
    smooth=True,
    region_bed=region_bed_file,
)

np.savez_compressed(sm.output["dp_mtx"], mat=dp_corrected)

# --- Window-level RDR ---
has_normal = "normal" in sample_types
logging.info(f"compute window RDRs, has_normal={has_normal}")

if has_normal:
    tumor_sidx = 1
    window_sizes = (bias_bed["END"] - bias_bed["START"]).to_numpy(dtype=np.float64)
    total_bases = np.nansum(dp_corrected * window_sizes[:, None], axis=0)
    library_correction = total_bases[0] / total_bases[1:]
    logging.info(f"library normalization factor: {library_correction}")

    normal_dp = dp_corrected[:, 0]
    valid = np.isfinite(normal_dp) & (normal_dp > 0)
    _n_valid = int(valid.sum())
    logging.info(f"normal: {_n_valid}/{n_windows} ({_n_valid / max(n_windows, 1) * 100:.1f}%) valid")

    with np.errstate(invalid="ignore", divide="ignore"):
        rdr_mat = dp_corrected[:, 1:] / normal_dp[:, None]
        rdr_mat *= library_correction[None, :]
    rdr_mat[~valid, :] = np.nan
else:
    tumor_sidx = 0
    rdr_mat = np.full_like(dp_corrected, np.nan)
    for i in range(nsamples):
        col = dp_corrected[:, i]
        valid_i = np.isfinite(col) & (col > 0)
        if not valid_i.any():
            logging.warning(
                f"  {rep_ids[i]}: no valid depth, skipping median-centering"
            )
            continue
        med = np.median(col[valid_i])
        logging.info(f"  median-centering {rep_ids[i]}: median={med:.4f}")
        with np.errstate(invalid="ignore", divide="ignore"):
            rdr_mat[valid_i, i] = col[valid_i] / med

tumor_rep_ids = rep_ids[tumor_sidx:]
n_nan_rdr = int(np.isnan(rdr_mat[:, 0]).sum())
logging.info(f"window RDR: {n_nan_rdr}/{n_windows} ({n_nan_rdr / max(n_windows, 1) * 100:.1f}%) NaN")

rdr_ylim = np.round(np.nanquantile(rdr_mat, 0.99)).astype(int) + 1
log_mad_and_plot(
    bias_bed,
    rdr_mat,
    tumor_rep_ids,
    genome_size,
    win_qc_dir,
    "rdr",
    "window",
    "RDR",
    rdr_ylim,
    smooth=True,
    region_bed=region_bed_file,
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

# --- Aggregate windows -> adaptive bins ---
logging.info("aggregating window data into adaptive bins")
bbs = pd.read_table(bb_file, sep="\t")
n_bb = len(bbs)
bb_dp = np.full((n_bb, nsamples), np.nan, dtype=np.float32)

bias_cols = ["GC"]
if "MAP" in bias_bed.columns:
    bias_cols.append("MAP")
if "REPLI" in bias_bed.columns:
    bias_cols.append("REPLI")
bias_arrays = {col: bias_bed[col].to_numpy(dtype=np.float64) for col in bias_cols}
bb_bias = {col: np.full(n_bb, np.nan, dtype=np.float64) for col in bias_cols}

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

    ov_win, ov_bin, ov_wt = compute_overlap_weights(
        win_start[win_idx_chr],
        win_end[win_idx_chr],
        bb_starts,
        bb_ends,
    )
    if len(ov_win) == 0:
        logging.info(f"  {chrom:<5s}: {n_bb_chr:>5d} bins, {0:>7d} windows")
        continue
    n_assigned = len(np.unique(ov_win))

    for s in range(nsamples):
        vals = dp_corrected[win_idx_chr, s][ov_win]
        bb_dp[bb_global_idx, s] = weighted_bincount_mean(vals, ov_bin, ov_wt, n_bb_chr)

    for col in bias_cols:
        vals = bias_arrays[col][win_idx_chr][ov_win]
        bb_bias[col][bb_global_idx] = weighted_bincount_mean(
            vals, ov_bin, ov_wt, n_bb_chr
        )

    logging.info(f"  {chrom:<5s}: {n_bb_chr:>5d} bins, {n_assigned:>7d} windows")

bb_corr = bbs[["#CHR", "START", "END"]].copy()
for col in bias_cols:
    bb_corr[col] = bb_bias[col]
bb_corr.to_csv(
    sm.output["corr_factors"], sep="\t", header=True, index=False, compression="gzip"
)
logging.info(f"wrote bb corr_factors: {bias_cols}, {n_bb} rows")

log_nan_summary("bb depth", bb_dp, rep_ids, n_bb)

bb_rd_ylim = max(np.nanquantile(bb_dp, 0.99), 1.0)
log_mad_and_plot(
    bbs,
    bb_dp,
    rep_ids,
    genome_size,
    bb_qc_dir,
    "depth",
    "bb",
    "RD",
    bb_rd_ylim,
    smooth=True,
    region_bed=region_bed_file,
)

for col in bias_cols:
    plot_1d_sample(
        bias_bed,
        bias_bed[col].to_numpy(),
        genome_size,
        os.path.join(win_qc_dir, f"{col}.pdf"),
        unit="window",
        val_type=col,
        smooth=True,
        region_bed=region_bed_file,
    )
    plot_1d_sample(
        bbs,
        bb_bias[col],
        genome_size,
        os.path.join(bb_qc_dir, f"{col}.pdf"),
        unit="bb",
        val_type=col,
        smooth=True,
        region_bed=region_bed_file,
    )

# --- Bin-level RDR ---
logging.info("computing bb RDR from aggregated depth")
if has_normal:
    normal_bb_dp = bb_dp[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        bb_rdr = bb_dp[:, tumor_sidx:] / normal_bb_dp[:, None]
        bb_rdr *= library_correction[None, :]
else:
    bb_rdr = np.full((n_bb, nsamples), np.nan, dtype=np.float32)
    for i in range(nsamples):
        col = bb_dp[:, i]
        valid_i = np.isfinite(col) & (col > 0)
        if valid_i.any():
            med = np.median(col[valid_i])
            logging.info(f"  bb median-centering {rep_ids[i]}: median={med:.4f}")
            with np.errstate(invalid="ignore", divide="ignore"):
                bb_rdr[valid_i, i] = col[valid_i] / med

n_nan_bb = int(np.isnan(bb_rdr[:, 0]).sum())
_n_filled = n_bb - n_nan_bb
logging.info(f"bb RDR: {_n_filled}/{n_bb} ({_n_filled / max(n_bb, 1) * 100:.1f}%) filled, {n_nan_bb} NaN")

np.savez_compressed(sm.output["rdr_mtx_bb"], mat=bb_rdr)
np.savez_compressed(sm.output["dp_mtx_bb"], mat=bb_dp)

log_mad_and_plot(
    bbs,
    bb_rdr,
    tumor_rep_ids,
    genome_size,
    bb_qc_dir,
    "rdr",
    "bb",
    "RDR",
    rdr_ylim,
    smooth=True,
    region_bed=region_bed_file,
)

logging.info("finished.")
