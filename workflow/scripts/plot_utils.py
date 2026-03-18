"""Plotting utilities for QC and diagnostic visualizations.
"""

import os
import logging

import numpy as np
import pandas as pd

from scipy.stats import gaussian_kde, pearsonr, spearmanr
from scipy.sparse import issparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import get_chr_sizes, read_region_file
from utils import adaptive_dot_size, stamp_path
from combine_counts_utils import compute_af_pseudobulk, compute_af_per_sample


# ---------------------------------------------------------------------------
# rd_correct_utils plots
# ---------------------------------------------------------------------------


def _plot_cov_panel(ax, covariate, reads, xlabel, show_ylabel, rep_id,
                    rmse=None, is_before=False, xlim=None, xticks=None):
    """KDE density scatter of readcov vs a covariate on a single axes."""
    valid = (reads > 0) & np.isfinite(covariate)
    if valid.sum() < 20:
        ax.set_visible(False)
        return
    x, y = covariate[valid], reads[valid]
    ylim = np.nanquantile(y, 0.99) * 1.1

    n_pts = len(x)
    rng = np.random.default_rng(0)
    if n_pts > 20000:
        idx = rng.choice(n_pts, size=20000, replace=False)
        kde = gaussian_kde(np.vstack([x[idx], y[idx]]))
    else:
        kde = gaussian_kde(np.vstack([x, y]))

    xlo = xlim[0] if xlim is not None else x.min()
    xhi = xlim[1] if xlim is not None else x.max()
    xgrid = np.linspace(xlo, xhi, 200)
    ygrid = np.linspace(0, ylim * 1.5, 200)
    xx, yy = np.meshgrid(xgrid, ygrid)
    positions = np.vstack([xx.ravel(), yy.ravel()])
    zz = kde(positions).reshape(xx.shape)

    ax.pcolormesh(xx, yy, zz, shading="gouraud", cmap="Blues", rasterized=True)
    ax.contour(
        xx, yy, zz, levels=6, colors="steelblue", linewidths=0.5, alpha=0.5
    )

    mad = np.median(np.abs(y - np.median(y)))
    r_pearson, _ = pearsonr(x, y)
    r_spearman, _ = spearmanr(x, y)
    metrics = f"MAD={mad:.2f}  r={r_pearson:.4f}  rho={r_spearman:.4f}"
    if rmse is not None and is_before:
        metrics = f"RMSE={rmse:.2f}  " + metrics
    logging.info(f"  {xlabel:<8s} {rep_id:<8s}: {metrics}")
    ax.set_xlabel(xlabel)
    if show_ylabel:
        ax.set_ylabel("Observed Readcov")
    ax.set_title(f"{rep_id}\n{metrics}", fontsize=8)
    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        ax.set_xlim(x.min(), x.max())
    if xticks is not None:
        ax.set_xticks(xticks)
    ax.set_ylim(0, ylim * 1.1)


def plot_gc_correction_pdf(gc, dp_before, dp_after, rep_ids, pdf, gc_rmse=None,
                           mappability=None, repliseq=None, title_prefix=""):
    """Two-page PDF: before/after correction KDE density plots.

    Each page has up to 3 rows (GC, MAP, RT) x nsamples columns.

    Parameters
    ----------
    gc_rmse : list of float or None
        Per-sample RMSE from the GC fit. Shown on the "Before" panel only.
    mappability : np.ndarray or None
        Per-window mappability values. If provided, a MAP row is added.
    repliseq : np.ndarray or None
        Per-window replication timing values. If provided, an RT row is added.
    title_prefix : str
        Optional prefix for page titles (e.g. ``"target — "``).
    """
    nsamples = len(rep_ids)
    panel_w = max(5, 5 * nsamples)

    # Build list of (row_label, covariate, xlabel)
    rows = [("GC", gc, "GC Content")]
    if mappability is not None:
        rows.append(("MAP", mappability, "Mappability"))
    if repliseq is not None:
        rows.append(("RT", repliseq, "Replication Timing"))
    nrows = len(rows)

    is_before = True
    for title, dp_mat in [
        ("Before Correction", dp_before),
        ("After Correction", dp_after),
    ]:
        logging.info(f"====={title}=====")
        fig, axes = plt.subplots(
            nrows, nsamples, figsize=(panel_w, 5 * nrows), squeeze=False,
        )
        for ri, (row_label, covariate, xlabel) in enumerate(rows):
            for si, rep_id in enumerate(rep_ids):
                rmse = gc_rmse[si] if (gc_rmse is not None and row_label == "GC") else None
                kw = {}
                if row_label == "MAP":
                    kw = {"xlim": (-0.2, 1.2), "xticks": np.arange(0, 1.1, 0.2)}
                _plot_cov_panel(
                    axes[ri, si],
                    covariate,
                    dp_mat[:, si],
                    xlabel,
                    show_ylabel=(si == 0),
                    rep_id=rep_id,
                    rmse=rmse,
                    is_before=is_before,
                    **kw,
                )
        fig.suptitle(f"{title_prefix}{title}", fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        is_before = False


# ---------------------------------------------------------------------------
# count_reads_utils plots
# ---------------------------------------------------------------------------


def plot_rd_gc(
    pos_df,
    mat,
    labels,
    genome_size,
    qc_dir,
    prefix,
    unit,
    val_type,
    ylim,
    gc_corr=None,
    gc_bin_median_std=None,
    run_id="",
    **plot_kwargs,
):
    """Generate a combined multi-sample 1-D chromosome scatter PDF.

    Produces a single ``{prefix}.pdf`` with one page per chromosome.  Each page
    contains vertically stacked subplots (one row per sample, shared x-axis).

    Parameters
    ----------
    gc_corr : dict[str, tuple[float, float]] | None
        {label: (pearson_r, spearman_r)} from ``compute_gc_rd_stats``.
    gc_bin_median_std : dict[str, float] | None
        {label: std of per-GC-bin median RD} from ``compute_gc_rd_stats``.
    """
    plot_file = stamp_path(os.path.join(qc_dir, f"{prefix}.pdf"), run_id)
    plot_1d_multi_sample(
        pos_df,
        mat,
        labels,
        genome_size,
        plot_file,
        unit=unit,
        val_type=val_type,
        max_ylim=ylim,
        **plot_kwargs,
    )


def plot_1d_multi_sample(
    pos_df: pd.DataFrame,
    mat: np.ndarray,
    labels: list,
    genome_size: str,
    out_file: str,
    unit="window",
    val_type="RD",
    s=4,
    dpi=72,
    alpha=0.6,
    min_ylim=0.0,
    max_ylim=None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
):
    """Multi-sample 1-D chromosome scatter plot: one page per chromosome,
    vertically stacked subplots (one row per sample, shared x-axis).

    Parameters
    ----------
    mat : np.ndarray
        (n_windows, n_samples) value matrix.
    labels : list[str]
        Sample labels, length == mat.shape[1].
    """
    n_samples = len(labels)
    logging.info(
        f"chromosome wide {unit}-level {val_type} multi-sample plot "
        f"({n_samples} samples), out_file={out_file}"
    )
    chrom_sizes = get_chr_sizes(genome_size)

    _region_by_chr = {}
    if region_bed is not None:
        _region_df = read_region_file(region_bed)
        for _ch, _grp in _region_df.groupby("#CHR", sort=False):
            _region_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    _blacklist_by_chr = {}
    if blacklist_bed is not None:
        _bl_df = read_region_file(blacklist_bed)
        for _ch, _grp in _bl_df.groupby("#CHR", sort=False):
            _blacklist_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    ch = pos_df["#CHR"].to_numpy()
    if "POS" in pos_df.columns:
        pos = pos_df["POS"].to_numpy()
    else:
        pos = ((pos_df["START"].to_numpy() + pos_df["END"].to_numpy()) // 2).astype(
            np.int64
        )
    num_bins = len(pos_df)
    change = np.flatnonzero(ch[1:] != ch[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, num_bins]
    chroms = ch[starts]

    pdf_fd = PdfPages(out_file)
    for chrom, lo, hi in zip(chroms, starts, ends):
        chr_end = chrom_sizes.get(chrom)
        if chr_end is None:
            logging.warning(f"{chrom}: not found in {genome_size}")
            continue

        x = pos[lo:hi]
        n_chr_pts = hi - lo
        s_chr = adaptive_dot_size(n_chr_pts, s_base=s)

        fig, axes = plt.subplots(
            nrows=n_samples,
            ncols=1,
            figsize=(40, 3 * n_samples),
            sharex=True,
            squeeze=False,
        )
        axes = axes[:, 0]

        for si, (ax, label) in enumerate(zip(axes, labels)):
            y = mat[lo:hi, si] if mat.ndim == 2 else mat[lo:hi]
            m = np.isfinite(y)

            if chrom in _region_by_chr:
                for reg_start, reg_end in _region_by_chr[chrom]:
                    ax.axvspan(
                        reg_start, reg_end, color="lightblue", alpha=0.15, zorder=0
                    )
            if chrom in _blacklist_by_chr:
                for bl_start, bl_end in _blacklist_by_chr[chrom]:
                    ax.axvspan(
                        bl_start, bl_end, color="lightcoral", alpha=0.2, zorder=0
                    )

            if m.any():
                ax.scatter(x[m], y[m], s=s_chr, alpha=alpha, rasterized=True)

            if val_type in ["AF", "BAF"]:
                ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
                ax.set_ylim(-0.05, 1.05)
            elif max_ylim is not None:
                ax.set_ylim(min_ylim, max_ylim)

            ax.set_ylabel(label, fontsize=10)
            ax.set_xlim(0, chr_end)
            if val_type not in ["AF", "BAF"]:
                ax.grid(axis="y", alpha=0.2)

            if si == 0:
                ax.set_title(f"{val_type} plot - {chrom}")

        axes[-1].set_xlabel(f"Position ({unit})")
        fig.tight_layout()
        pdf_fd.savefig(fig, dpi=dpi)
        plt.close(fig)
    pdf_fd.close()


def plot_1d_sample(
    pos_df: pd.DataFrame,
    val: np.ndarray,
    genome_size: str,
    out_file: str,
    unit="SNP",
    val_type="BAF",
    s=4,
    dpi=72,
    alpha=0.6,
    figsize=(40, 3),
    min_ylim=0.0,
    max_ylim=None,
    mask: np.ndarray | None = None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
):
    """
    plot any features like AF, RDR, etc., in 1d chromosome scatter plot.
    """
    logging.info(f"chromosome wide {unit}-level {val_type} plot, out_file={out_file}")
    chrom_sizes = get_chr_sizes(genome_size)

    _region_by_chr = {}
    if region_bed is not None:
        _region_df = read_region_file(region_bed)
        for _ch, _grp in _region_df.groupby("#CHR", sort=False):
            _region_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    _blacklist_by_chr = {}
    if blacklist_bed is not None:
        _bl_df = read_region_file(blacklist_bed)
        for _ch, _grp in _bl_df.groupby("#CHR", sort=False):
            _blacklist_by_chr[_ch] = list(zip(_grp["START"], _grp["END"]))

    ch = pos_df["#CHR"].to_numpy()
    if "POS" in pos_df.columns:
        pos = pos_df["POS"].to_numpy()
    else:
        pos = ((pos_df["START"].to_numpy() + pos_df["END"].to_numpy()) // 2).astype(
            np.int64
        )
    num_bins = len(pos_df)
    change = np.flatnonzero(ch[1:] != ch[:-1]) + 1
    starts = np.r_[0, change]
    ends = np.r_[change, num_bins]
    chroms = ch[starts]

    pdf_fd = PdfPages(out_file)
    for chrom, lo, hi in zip(chroms, starts, ends):
        chr_end = chrom_sizes.get(chrom)
        if chr_end is None:
            logging.warning(f"{chrom}: not found in {genome_size}")
            continue

        x = pos[lo:hi]
        y = val[lo:hi]
        m = np.isfinite(y)
        mask_chr = mask[lo:hi] if mask is not None else None

        n_plot = int(m.sum())
        s_chr = adaptive_dot_size(n_plot, s_base=s)
        if n_plot == 0:
            logging.warning(f"{chrom}: all {val_type} values are non-finite")
            continue
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        if chrom in _region_by_chr:
            for reg_start, reg_end in _region_by_chr[chrom]:
                ax.axvspan(reg_start, reg_end, color="lightblue", alpha=0.15, zorder=0)
        if chrom in _blacklist_by_chr:
            for bl_start, bl_end in _blacklist_by_chr[chrom]:
                ax.axvspan(bl_start, bl_end, color="lightcoral", alpha=0.2, zorder=0)
        if mask_chr is not None:
            kept = m & mask_chr
            filt = m & ~mask_chr
            if filt.any():
                ax.scatter(
                    x[filt],
                    y[filt],
                    s=s_chr,
                    alpha=0.8,
                    color="red",
                    rasterized=True,
                    label=f"filtered ({filt.sum()})",
                )
            if kept.any():
                ax.scatter(
                    x[kept],
                    y[kept],
                    s=s_chr,
                    alpha=0.8,
                    color="blue",
                    rasterized=True,
                    label=f"kept ({kept.sum()})",
                )
            ax.legend(loc="upper right", fontsize=8, markerscale=2)
        else:
            ax.scatter(x[m], y[m], s=s_chr, alpha=alpha, rasterized=True)
        if val_type in ["AF", "BAF"]:
            ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
            ax.set_ylim(-0.05, 1.05)
        elif max_ylim is not None:
            ax.set_ylim(min_ylim, max_ylim)

        ax.set_xlabel(val_type)
        ax.set_xlim(0, chr_end)
        if val_type not in ["AF", "BAF"]:
            ax.grid(axis="y", alpha=0.2)
        ax.set_title(f"{val_type} plot - {chrom}")
        pdf_fd.savefig(fig, dpi=dpi)
        plt.close(fig)
    pdf_fd.close()
    return


# ---------------------------------------------------------------------------
# combine_counts_utils plots
# ---------------------------------------------------------------------------


def plot_snp_depth_histogram(
    tot_mtx,
    rep_ids,
    qc_dir,
    run_id,
    ref_mtx=None,
    is_bulk=True,
):
    """Plot per-sample histograms of total allele depth and ref-AF at SNP positions.

    The caller is responsible for pre-slicing *tot_mtx* (and *ref_mtx*) to the
    desired SNP rows before calling this function.

    Parameters
    ----------
    tot_mtx : ndarray or sparse
        Total depth matrix (SNPs x samples/cells), already filtered to the
        SNPs of interest.
    rep_ids : list[str]
        Sample / replicate identifiers.
    qc_dir : str
        Output directory for the PDF.
    run_id : str
        Run identifier used by :func:`stamp_path`.
    ref_mtx : ndarray or sparse, optional
        Ref-allele count matrix (same shape as *tot_mtx*). When provided, a
        second column of ref-AF histograms is added to the figure.
    is_bulk : bool
        If True, one row per sample. If False, aggregate across cells into a
        single pseudobulk row.
    """
    logging.info("QC analysis - plot SNP depth histogram")

    has_af = ref_mtx is not None
    ncols = 2 if has_af else 1

    if is_bulk:
        nrows = len(rep_ids)
        row_labels = list(rep_ids)
    else:
        nrows = 1
        row_labels = ["pseudobulk"]

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(5 * ncols, 3 * nrows), squeeze=False
    )

    def _stats_title(name, depth_vals):
        if len(depth_vals) == 0:
            return f"{name} (n=0)"
        return (
            f"{name}\n"
            f"mean={depth_vals.mean():.1f}, median={np.median(depth_vals):.1f}, "
            f"min={depth_vals.min():.0f}, max={depth_vals.max():.0f}"
        )

    def _get_col(mat, col_idx):
        """Extract a single column as a 1-D numpy array."""
        if issparse(mat):
            return np.asarray(mat[:, col_idx].toarray()).ravel()
        return np.asarray(mat[:, col_idx]).ravel()

    def _pseudobulk_sum(mat):
        """Sum across all columns (cells) to produce a 1-D pseudobulk array."""
        if issparse(mat):
            return np.asarray(mat.sum(axis=1)).ravel()
        return np.asarray(mat.sum(axis=1)).ravel()

    for ri, label in enumerate(row_labels):
        # --- depth histogram ---
        ax_depth = axes[ri, 0]
        if is_bulk:
            depth = _get_col(tot_mtx, ri)
        else:
            depth = _pseudobulk_sum(tot_mtx)
        if len(depth) > 0:
            clip_threshold = np.percentile(depth, 99)
            ax_depth.hist(depth[depth <= clip_threshold], bins=50, alpha=0.7)
        ax_depth.set_title(_stats_title(label, depth), fontsize=9)
        ax_depth.set_xlabel("Total allele depth")
        ax_depth.set_ylabel("# SNPs")

        # --- ref-AF histogram ---
        if has_af:
            ax_af = axes[ri, 1]
            if is_bulk:
                total_depth = _get_col(tot_mtx, ri).astype(np.float64)
                ref_depth = _get_col(ref_mtx, ri).astype(np.float64)
            else:
                total_depth = _pseudobulk_sum(tot_mtx).astype(np.float64)
                ref_depth = _pseudobulk_sum(ref_mtx).astype(np.float64)
            covered_mask = total_depth > 0
            ref_af = np.full(len(total_depth), np.nan)
            ref_af[covered_mask] = ref_depth[covered_mask] / total_depth[covered_mask]
            ref_af_valid = ref_af[~np.isnan(ref_af)]
            if len(ref_af_valid) > 0:
                ax_af.hist(ref_af_valid, bins=50, range=(0, 1), alpha=0.7)
            ax_af.axvline(0.5, color="red", linestyle="--", linewidth=0.8)
            ax_af.set_xlim(0, 1)
            ax_af.set_title(f"{label} ref-AF", fontsize=9)
            ax_af.set_xlabel("Ref allele frequency")
            ax_af.set_ylabel("# SNPs")

    fig.tight_layout()
    out_path = stamp_path(os.path.join(qc_dir, "snp_depth_hist.pdf"), run_id)
    fig.savefig(out_path)
    plt.close(fig)
    logging.info(f"saved SNP depth histogram to {out_path}")


def plot_allele_freqs(
    pos_df,
    rep_ids,
    tot_mtx,
    b_mtx,
    genome_size,
    plot_dir,
    apply_pseudobulk=False,
    allele="ref",
    unit="SNP",
    suffix="",
    snp_mask=None,
    region_bed=None,
    blacklist_bed=None,
    run_id="",
):
    """Generate genome-wide allele-frequency scatter plots.

    Produces either one plot per sample or a single pseudobulk plot, saved as
    multi-page PDFs (one page per chromosome).

    Parameters
    ----------
    pos_df : pd.DataFrame
        SNP/bin position DataFrame with ``#CHR`` and ``POS`` (or ``START``/``END``).
    rep_ids : list[str]
        Replicate identifiers.
    tot_mtx, b_mtx : sparse or ndarray
        Total depth and B-allele count matrices.
    genome_size : str
        Path to chromosome sizes file.
    plot_dir : str
        Output directory for PDF plots.
    apply_pseudobulk : bool
        If True, plot a single pseudobulk AF; otherwise plot per-sample.
    allele : str
        Allele label for filenames (e.g., ``"ref"``, ``"B"``).
    unit : str
        Feature unit label (e.g., ``"SNP"``, ``"bin"``).
    suffix : str
        Optional filename suffix.
    region_bed : str or None
        Path to whitelist region BED file for background shading.
    blacklist_bed : str or None
        Path to blacklist BED file for background shading.
    """
    logging.info(
        f"QC analysis - plot {allele}-{unit} allele frequency, {unit}={allele}, apply_pseudobulk={apply_pseudobulk}"
    )
    if apply_pseudobulk:
        af = compute_af_pseudobulk(tot_mtx, b_mtx)
        plot_file = stamp_path(os.path.join(plot_dir, f"af_{allele}_{unit}.pseudobulk{suffix}.pdf"), run_id)
        plot_1d_sample(
            pos_df,
            af,
            genome_size,
            plot_file,
            unit=unit,
            val_type="AF",
            mask=snp_mask,
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
        )
    else:
        _tot_mtx = tot_mtx.tocsc() if issparse(tot_mtx) else tot_mtx
        _b_mtx = b_mtx.tocsc() if issparse(b_mtx) else b_mtx
        af_mat = np.column_stack(
            [compute_af_per_sample(_tot_mtx, _b_mtx, i) for i in range(len(rep_ids))]
        )
        plot_file = stamp_path(os.path.join(plot_dir, f"af_{allele}_{unit}{suffix}.pdf"), run_id)
        plot_1d_multi_sample(
            pos_df,
            af_mat,
            list(rep_ids),
            genome_size,
            plot_file,
            unit=unit,
            val_type="AF",
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
        )
    return
