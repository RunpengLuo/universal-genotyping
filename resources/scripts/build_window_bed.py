#!/usr/bin/env python3
"""Build a per-window bias correction BED with region/blacklist/WES filtering.

Generates windows (WGS fixed-size or WES target+antitarget), filters by
accessible regions and blacklist, then computes per-window GC content with
optional mappability and replication timing columns.

WGS mode (default):
  1. Subtract --blacklist_bed from --region_bed (if provided)
  2. Tile directly within each remaining interval (every window fully contained)
  3. Merge undersized trailing windows into the previous window

WES mode (--wes_targets_bed):
  1. Split vendor exon targets into ~target_avg_size pieces
  2. Generate antitarget bins in accessible non-target regions
  3. Combine target + antitarget, mark is_target column

Both modes then subtract blacklist, assign region_id, and compute GC/MAP/REPLI.

Output is a gzipped TSV:
  #CHR  START  END  region_id  GC  [MAP]  [REPLI]  [is_target]

Dependencies:
  - pybedtools, pyranges, pandas, numpy
  - bigWigToBedGraph, liftOver (UCSC tools) -- only if --repliseq
"""

import argparse
import glob
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd
import pyranges as pr
from pybedtools import BedTool

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHROM_ORDER = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]

UCSC_BASE = (
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq"
)
WAVE_SIGNAL_FILES = [
    "wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig",
    "wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig",
]

CHAIN_URL = (
    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_chrom_sizes(genome_size_file):
    """Read chromosome sizes file into a DataFrame."""
    return pd.read_csv(
        genome_size_file,
        sep="\t",
        header=None,
        names=["chrom", "size"],
        dtype={"chrom": str, "size": np.int64},
    )


def _df_to_pyranges(df):
    """Convert a DataFrame with #CHR/START/END to PyRanges with _idx.

    Uses a DataFrame constructor so that _idx values are preserved through
    pyranges' internal chromosome reordering (insert() renumbers them).
    """
    pr_df = pd.DataFrame({
        "Chromosome": df["#CHR"].values,
        "Start": df["START"].values,
        "End": df["END"].values,
        "_idx": np.arange(len(df)),
    })
    return pr.PyRanges(pr_df)


def _subdivide_intervals(df, chrom_col, start_col, end_col, avg_size, min_size=0):
    """Subdivide intervals into ~avg_size pieces.

    Returns list of (chrom, start, end) tuples.
    """
    rows = []
    for _, r in df.iterrows():
        span = r[end_col] - r[start_col]
        nbins = max(1, round(span / avg_size))
        bin_size = span / nbins
        for i in range(nbins):
            s = r[start_col] + int(round(i * bin_size))
            e = r[start_col] + int(round((i + 1) * bin_size))
            if (e - s) >= min_size:
                rows.append((r[chrom_col], s, e))
    return rows


def _sort_windows(windows):
    """Sort windows by natural chromosome order, then by start position."""
    chrom_rank = {c: i for i, c in enumerate(CHROM_ORDER)}
    windows = windows[windows["#CHR"].isin(chrom_rank)].copy()
    windows["_chrom_rank"] = windows["#CHR"].map(chrom_rank)
    windows = windows.sort_values(["_chrom_rank", "START"]).reset_index(drop=True)
    return windows.drop(columns=["_chrom_rank"])


# ---------------------------------------------------------------------------
# Reference version helpers
# ---------------------------------------------------------------------------


def get_standard_chroms(reference_version, genome_size_file):
    """Derive valid standard chromosomes from reference version and genome size file."""
    all_chroms = set(_read_chrom_sizes(genome_size_file)["chrom"])
    candidates = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
    if reference_version not in ("hg19", "hg38", "chm13v2"):
        candidates |= {str(c) for c in list(range(1, 23)) + ["X", "Y"]}
    return candidates & all_chroms


# ---------------------------------------------------------------------------
# Window generation (WGS)
# ---------------------------------------------------------------------------


def _tile_region(chrom, start, end, window_size):
    """Tile a single region into fixed-size windows, merging last short bin.

    If the last window is shorter than window_size // 2 and there is a previous
    window, it is merged into the previous window (extending its END).  Regions
    smaller than window_size // 2 still emit one window.
    """
    rows = []
    pos = start
    while pos < end:
        w_end = min(pos + window_size, end)
        rows.append([chrom, pos, w_end])
        pos = w_end
    # Merge last bin into previous if undersized
    if len(rows) > 1 and (rows[-1][2] - rows[-1][1]) < window_size // 2:
        rows[-2][2] = rows[-1][2]
        rows.pop()
    return rows


def generate_wgs_windows(genome_size_file, window_size, standard_chroms,
                         region_bed, blacklist_bed):
    """Tile windows directly within accessible regions.

    Reads the region BED, subtracts the blacklist, then tiles each remaining
    interval into fixed-size windows (merging undersized trailing windows).
    Every output window is fully contained in a blacklist-subtracted region.
    """
    regions = pr.read_bed(region_bed)
    if blacklist_bed:
        regions = regions.subtract(pr.read_bed(blacklist_bed))

    reg_df = regions.df
    reg_df = reg_df[reg_df["Chromosome"].isin(standard_chroms)].reset_index(drop=True)

    rows = []
    for _, r in reg_df.iterrows():
        rows.extend(_tile_region(r["Chromosome"], r["Start"], r["End"], window_size))

    windows = pd.DataFrame(rows, columns=["#CHR", "START", "END"])
    print(f"  {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes")
    return windows


# ---------------------------------------------------------------------------
# Window generation (WES)
# ---------------------------------------------------------------------------


def generate_wes_windows(wes_targets_bed, region_bed, blacklist_bed,
                         target_avg_size, antitarget_avg_size, pad_size,
                         standard_chroms):
    """Generate target + antitarget windows for WES mode."""
    # --- Targets: subdivide vendor exon capture intervals ---
    targets = _read_and_subdivide_targets(wes_targets_bed, target_avg_size)
    targets["is_target"] = 1

    # --- Antitargets: off-target bins in accessible non-exon regions ---
    antitargets = _generate_antitargets(
        targets, region_bed, blacklist_bed, antitarget_avg_size, pad_size,
    )
    antitargets["is_target"] = 0

    windows = pd.concat([targets, antitargets], ignore_index=True)
    windows = windows[windows["#CHR"].isin(standard_chroms)].reset_index(drop=True)

    n_tgt = int(windows["is_target"].sum())
    print(f"  {len(windows)} windows ({n_tgt} target, {len(windows) - n_tgt} antitarget)")
    return windows


def _dedup_targets(df):
    """Deduplicate overlapping targets CNVkit-style.

    1. Exact duplicates: keep first.
    2. Partial overlap [a, b]: truncate b's start to a's end.
    3. Target a fully inside target b: skip a.

    Assumes df is sorted by (#CHR, START).  Ties in START are broken by
    descending END so the largest interval is processed first.
    """
    df = df.sort_values(
        ["#CHR", "START", "END"], ascending=[True, True, False],
    ).reset_index(drop=True)

    rows = []
    prev_chrom = None
    prev_end = 0
    n_dup, n_contained, n_truncated = 0, 0, 0

    for _, r in df.iterrows():
        chrom, start, end = r["#CHR"], r["START"], r["END"]

        if chrom != prev_chrom:
            rows.append((chrom, start, end))
            prev_chrom, prev_end = chrom, end
            continue

        if start == rows[-1][1] and end == rows[-1][2]:
            n_dup += 1
            continue

        if end <= prev_end:
            n_contained += 1
            continue

        if start < prev_end:
            n_truncated += 1
            start = prev_end

        if start < end:
            rows.append((chrom, start, end))
            prev_end = end

    print(
        f"  dedup: {n_dup} exact duplicates, {n_contained} contained, "
        f"{n_truncated} truncated"
    )
    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


def _read_and_subdivide_targets(wes_targets_bed, target_avg_size):
    """Read vendor exon targets, deduplicate CNVkit-style, subdivide."""
    df = pd.read_csv(
        wes_targets_bed, sep="\t", header=None,
        usecols=[0, 1, 2], names=["#CHR", "START", "END"], comment="#",
    )
    df = df[df["END"] > df["START"]].copy()
    n_raw = len(df)
    df = _dedup_targets(df)
    print(f"  {n_raw} raw targets, {len(df)} after dedup")
    rows = _subdivide_intervals(df, "#CHR", "START", "END", target_avg_size, min_size=1)
    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


def _generate_antitargets(targets, region_bed, blacklist_bed,
                          antitarget_avg_size, pad_size):
    """Generate off-target bins in accessible regions not covered by exon targets."""
    access = pr.read_bed(region_bed)
    if blacklist_bed:
        access = access.subtract(pr.read_bed(blacklist_bed))

    tgt_pr = pr.PyRanges(
        chromosomes=targets["#CHR"], starts=targets["START"], ends=targets["END"],
    )
    padded = tgt_pr.copy()
    padded.Start = (padded.Start - pad_size).clip(lower=0)
    padded.End = padded.End + pad_size

    background = access.subtract(padded)
    bg_df = background.df
    bg_df = bg_df[bg_df["End"] > bg_df["Start"]].copy()

    min_bin_size = antitarget_avg_size // 16
    rows = _subdivide_intervals(
        bg_df, "Chromosome", "Start", "End", antitarget_avg_size, min_size=min_bin_size,
    )
    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


# ---------------------------------------------------------------------------
# Blacklist subtraction
# ---------------------------------------------------------------------------


def subtract_blacklist(windows, blacklist_bed_file):
    """Remove windows overlapping blacklisted regions."""
    bl_idx = (
        _df_to_pyranges(windows)
        .overlap(pr.read_bed(blacklist_bed_file))
        .df["_idx"]
        .unique()
    )
    keep_mask = ~np.isin(np.arange(len(windows)), bl_idx)
    return windows.loc[keep_mask].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Region ID assignment
# ---------------------------------------------------------------------------


def assign_region_id(windows, region_bed_file):
    """Assign region_id to each window by mapping midpoints to region intervals.

    Returns region_id as strings in ``CHR:START-END`` format, or NaN.
    """
    regions = pr.read_bed(region_bed_file).df
    regions = regions.rename(
        columns={"Chromosome": "#CHR", "Start": "START", "End": "END"}
    )
    regions["region_id"] = (
        regions["#CHR"].astype(str)
        + ":"
        + regions["START"].astype(str)
        + "-"
        + regions["END"].astype(str)
    )

    mids = (windows["START"] + windows["END"]) // 2
    result = pd.Series(pd.NA, index=windows.index, dtype="object")

    for chrom in regions["#CHR"].unique():
        win_mask = windows["#CHR"] == chrom
        if not win_mask.any():
            continue
        reg_ch = regions.loc[regions["#CHR"] == chrom].sort_values("START")
        starts = reg_ch["START"].to_numpy()
        ends = reg_ch["END"].to_numpy()
        ids = reg_ch["region_id"].to_numpy()
        positions = mids[win_mask].to_numpy()

        idx = np.searchsorted(starts, positions, side="right") - 1
        safe_idx = idx.clip(min=0)
        valid = (idx >= 0) & (positions < ends[safe_idx])
        win_indices = windows.index[win_mask]
        result.loc[win_indices[valid]] = ids[idx[valid]]

    return result


# ---------------------------------------------------------------------------
# Per-window covariates
# ---------------------------------------------------------------------------


def compute_gc(windows_df, reference):
    """Compute per-window GC fraction via pybedtools nucleotide_content."""
    bt = BedTool.from_dataframe(windows_df[["#CHR", "START", "END"]])
    nuc = bt.nucleotide_content(fi=reference).to_dataframe(disable_auto_names=True)
    windows_df = windows_df.copy()
    windows_df["GC"] = nuc["5_pct_gc"].values
    return windows_df


def compute_mappability(windows_df, mappability_file, genome_size_file):
    """Compute per-window mean mappability via pybedtools map."""
    # Add _idx so we can restore original row order after bedtools sort
    win_bed = windows_df[["#CHR", "START", "END"]].copy()
    win_bed["_idx"] = np.arange(len(win_bed))
    bt = BedTool.from_dataframe(win_bed).sort(g=genome_size_file)
    map_bt = BedTool(mappability_file).sort(g=genome_size_file)
    map_cov = bt.map(b=map_bt, c=4, o="mean", g=genome_size_file).to_dataframe(
        disable_auto_names=True
    )
    map_cov.columns = ["#CHR", "START", "END", "_idx", "MAP"]
    map_cov["_idx"] = map_cov["_idx"].astype(int)
    map_cov = map_cov.sort_values("_idx")
    return pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(0.0).clip(0.0, 1.0).values


# ---------------------------------------------------------------------------
# Replication timing (Repli-seq)
# ---------------------------------------------------------------------------


def _download(url, dest):
    """Download a file with wget."""
    subprocess.check_call(["wget", "-q", "-O", dest, url])


def _bigwig_to_bedgraph(bigwig, bedgraph):
    """Convert bigWig to bedGraph using UCSC bigWigToBedGraph."""
    subprocess.check_call(["bigWigToBedGraph", bigwig, bedgraph])


def _liftover(input_bed, chain, output_bed, unmapped):
    """Run UCSC liftOver."""
    subprocess.check_call(["liftOver", input_bed, chain, output_bed, unmapped])


def _prepare_lifted_bedgraphs(bigwig_dir, work_dir, chain):
    """Download bigWigs if needed, convert to bedGraph, liftOver to hg38.

    Returns list of lifted bedGraph file paths.
    """
    os.makedirs(bigwig_dir, exist_ok=True)
    n_files = len(WAVE_SIGNAL_FILES)
    for i, fname in enumerate(WAVE_SIGNAL_FILES, 1):
        dest = os.path.join(bigwig_dir, fname)
        if not os.path.isfile(dest):
            print(f"\r  downloading bigWigs [{i}/{n_files}]", end="", flush=True)
            _download(f"{UCSC_BASE}/{fname}", dest)
    print()

    bigwig_files = sorted(glob.glob(os.path.join(bigwig_dir, "*WaveSignal*.bigWig")))
    if not bigwig_files:
        sys.exit(f"Error: no WaveSignal bigWig files found in {bigwig_dir}")

    lifted_dir = os.path.join(work_dir, "lifted")
    os.makedirs(lifted_dir, exist_ok=True)

    lifted_files = []
    for i, bw in enumerate(bigwig_files, 1):
        name = os.path.basename(bw).replace(".bigWig", "")
        bg_file = os.path.join(lifted_dir, f"{name}.hg19.bedGraph")
        lifted_file = os.path.join(lifted_dir, f"{name}.hg38.bedGraph")
        unmapped_file = os.path.join(lifted_dir, f"{name}.unmapped")

        if not os.path.isfile(lifted_file):
            print(f"\r  liftOver [{i}/{len(bigwig_files)}]", end="", flush=True)
            _bigwig_to_bedgraph(bw, bg_file)
            _liftover(bg_file, chain, lifted_file, unmapped_file)
            os.remove(bg_file)
        lifted_files.append(lifted_file)
    print()
    return lifted_files


def _bin_bedgraph_signals(windows_df, lifted_files, standard_chroms):
    """Bin lifted bedGraph signals into windows. Returns per-window mean signal."""
    n_win = len(windows_df)
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int32)

    for i, lf in enumerate(lifted_files, 1):
        print(f"\r  binning repli-seq [{i}/{len(lifted_files)}]", end="", flush=True)
        bg = pd.read_csv(
            lf, sep="\t", header=None,
            names=["chrom", "start", "end", "signal"],
            dtype={"chrom": str, "start": np.int64, "end": np.int64, "signal": np.float64},
        )
        bg = bg[bg["chrom"].isin(standard_chroms)].reset_index(drop=True)

        for chrom, grp in bg.groupby("chrom", sort=False):
            win_mask = windows_df["#CHR"] == chrom
            win_idx = np.where(win_mask)[0]
            if len(win_idx) == 0:
                continue

            win_starts = windows_df["START"].to_numpy()[win_idx]
            sort_order = np.argsort(win_starts)
            win_starts = win_starts[sort_order]
            win_idx = win_idx[sort_order]

            bg_mids = ((grp["start"] + grp["end"]) // 2).to_numpy()
            bin_idx = np.searchsorted(win_starts, bg_mids, side="right") - 1
            valid = (bin_idx >= 0) & (bin_idx < len(win_idx))
            global_idx = win_idx[bin_idx[valid]]
            signal_sum[global_idx] += grp["signal"].to_numpy()[valid]
            signal_count[global_idx] += 1
    print()

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_signal = np.where(signal_count > 0, signal_sum / signal_count, np.nan)

    n_covered = int((signal_count > 0).sum())
    print(f"  {n_covered}/{n_win} windows have replication timing data")
    return np.round(mean_signal, 6)


def compute_repliseq(windows_df, standard_chroms, args):
    """Compute per-window mean replication timing from ENCODE Repli-seq bigWigs.

    Handles work_dir/chain setup, then delegates to helpers.
    """
    if args.work_dir:
        work_dir = args.work_dir
        os.makedirs(work_dir, exist_ok=True)
        cleanup = False
    else:
        work_dir = tempfile.mkdtemp(prefix="repliseq_")
        cleanup = True

    chain = args.chain
    if chain is None:
        chain = os.path.join(work_dir, "hg19ToHg38.over.chain.gz")
        if not os.path.isfile(chain):
            print("  Downloading chain file ...")
            _download(CHAIN_URL, chain)
    if not os.path.isfile(chain):
        sys.exit(f"Error: chain file not found: {chain}")

    bigwig_dir = args.bigwig_dir or os.path.join(work_dir, "bigwigs")
    lifted_files = _prepare_lifted_bedgraphs(bigwig_dir, work_dir, chain)
    result = _bin_bedgraph_signals(windows_df, lifted_files, standard_chroms)

    if cleanup:
        import shutil
        shutil.rmtree(work_dir)

    return result


# ---------------------------------------------------------------------------
# Summary reporting
# ---------------------------------------------------------------------------


def print_summary(windows):
    """Print per-window and per-chromosome summary statistics."""
    sizes = windows["END"] - windows["START"]
    has_target = "is_target" in windows.columns

    print()
    print("  Summary:")
    print(f"    {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes")
    print(f"    {windows['region_id'].nunique()} distinct regions")
    print(
        f"    window size: min={sizes.min()}, median={int(sizes.median())}, "
        f"max={sizes.max()}, mean={sizes.mean():.0f}"
    )
    total_bp = int(sizes.sum())
    print(f"    total coverage: {total_bp:,} bp ({total_bp / 1e9:.2f} Gb)")

    if has_target:
        tgt_sizes = sizes[windows["is_target"].astype(bool)]
        anti_sizes = sizes[~windows["is_target"].astype(bool)]
        print(
            f"    target size:     min={tgt_sizes.min()}, median={int(tgt_sizes.median())}, "
            f"max={tgt_sizes.max()}"
        )
        print(
            f"    antitarget size: min={anti_sizes.min()}, median={int(anti_sizes.median())}, "
            f"max={anti_sizes.max()}"
        )

    # Per-chromosome table
    header = f"    {'chrom':<8} {'windows':>8} {'regions':>8} {'coverage_Mb':>12} {'median_size':>12}"
    if has_target:
        header += f" {'target':>8} {'antitarget':>10}"
    print(header)

    for chrom in CHROM_ORDER:
        mask = windows["#CHR"] == chrom
        if not mask.any():
            continue
        ch_sizes = sizes[mask]
        ch_regions = windows.loc[mask, "region_id"].nunique()
        row = (
            f"    {chrom:<8} {int(mask.sum()):>8} {ch_regions:>8} "
            f"{ch_sizes.sum() / 1e6:>12.2f} {int(ch_sizes.median()):>12}"
        )
        if has_target:
            n_tgt = int(windows.loc[mask, "is_target"].astype(bool).sum())
            n_anti = int(mask.sum()) - n_tgt
            row += f" {n_tgt:>8} {n_anti:>10}"
        print(row)
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a per-window bias correction BED. "
            "Generates windows with region/blacklist filtering (WGS) or "
            "target+antitarget splitting (WES), then computes per-window "
            "GC content with optional mappability and replication timing."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # WGS: GC-only (1 kb windows)\n"
            "  python build_window_bed.py \\\n"
            "    --reference_version hg38 \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --region_bed /path/to/regions.bed \\\n"
            "    --out_file gc.1kbp.hg38.bed.gz\n\n"
            "  # WES: target+antitarget with GC+MAP\n"
            "  python build_window_bed.py \\\n"
            "    --reference_version hg38 \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --region_bed /path/to/regions.bed \\\n"
            "    --wes_targets_bed /path/to/targets.bed \\\n"
            "    --mappability_bed /path/to/mappability.bed \\\n"
            "    --out_file gc_map.wes.bed.gz"
        ),
    )
    # Required
    parser.add_argument("--reference_version", required=True,
                        help="Reference genome version (hg19, hg38, chm13v2).")
    parser.add_argument("--reference", required=True,
                        help="Path to reference FASTA (must have .fai index).")
    parser.add_argument("--genome_size", required=True,
                        help="Path to chromosome sizes file (tab-separated: chrom, size).")
    parser.add_argument("--region_bed", required=True,
                        help="Accessible regions BED -- keep only windows overlapping these regions.")
    parser.add_argument("--out_file", required=True,
                        help="Output gzipped TSV path (e.g. gc_map_repli.1kbp.hg38.bed.gz).")
    # Optional filtering
    parser.add_argument("--blacklist_bed", default=None,
                        help="Optional blacklist BED -- subtract these regions.")
    # WGS options
    parser.add_argument("--window", type=int, default=1000,
                        help="Window size in bp (WGS only, default: 1000).")
    # WES options
    parser.add_argument("--wes_targets_bed", default=None,
                        help="Vendor exon capture targets BED. Enables WES mode.")
    parser.add_argument("--target_avg_size", type=int, default=200,
                        help="Average target bin size (WES only, default: 200).")
    parser.add_argument("--antitarget_avg_size", type=int, default=150000,
                        help="Average off-target bin size (WES only, default: 150000).")
    parser.add_argument("--pad_size", type=int, default=500,
                        help="Padding around targets for antitarget generation (WES only, default: 500).")
    # Covariate options
    parser.add_argument("--mappability_bed", default=None,
                        help="Optional BED-format mappability track (4th column = score).")
    parser.add_argument("--repliseq", action="store_true",
                        help="Enable replication timing from ENCODE Repli-seq bigWigs.")
    parser.add_argument("--chain", default=None,
                        help="hg19-to-hg38 liftOver chain file (downloaded if omitted with --repliseq).")
    parser.add_argument("--bigwig_dir", default=None,
                        help="Directory for pre-downloaded WaveSignal bigWig files.")
    parser.add_argument("--work_dir", default=None,
                        help="Working directory for intermediate files (default: temp dir).")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    # --- Validate inputs ---
    for name, path in [
        ("reference", args.reference),
        ("genome_size", args.genome_size),
        ("region_bed", args.region_bed),
        ("blacklist_bed", args.blacklist_bed),
        ("mappability_bed", args.mappability_bed),
        ("wes_targets_bed", args.wes_targets_bed),
    ]:
        if path and not os.path.isfile(path):
            sys.exit(f"Error: {name} not found: {path}")

    wes_mode = args.wes_targets_bed is not None
    t0 = time.time()

    print("=== build_window_bed.py ===")
    print(f"  Ref version: {args.reference_version}")
    print(f"  Reference:   {args.reference}")
    print(f"  Genome size: {args.genome_size}")
    print(f"  Region BED:  {args.region_bed}")
    print(f"  Blacklist:   {args.blacklist_bed or 'none'}")
    print(f"  Mode:        {'WES' if wes_mode else 'WGS'}")
    if wes_mode:
        print(f"  WES targets: {args.wes_targets_bed}")
        print(f"  Target avg:  {args.target_avg_size} bp")
        print(f"  Antitarget:  {args.antitarget_avg_size} bp")
        print(f"  Pad size:    {args.pad_size} bp")
    else:
        print(f"  Window:      {args.window} bp")
    print(f"  Mappability: {args.mappability_bed or 'none'}")
    print(f"  Repli-seq:   {'yes' if args.repliseq else 'no'}")
    print(f"  Output:      {args.out_file}")
    print()

    standard_chroms = get_standard_chroms(args.reference_version, args.genome_size)
    print(f"  {len(standard_chroms)} standard chromosomes")

    # --- Step 1: Generate windows ---
    if wes_mode:
        print("[1/6] WES: generating target + antitarget windows ...")
        windows = generate_wes_windows(
            args.wes_targets_bed, args.region_bed, args.blacklist_bed,
            args.target_avg_size, args.antitarget_avg_size, args.pad_size,
            standard_chroms,
        )
    else:
        print("[1/6] WGS: tiling windows within regions ...")
        windows = generate_wgs_windows(
            args.genome_size, args.window, standard_chroms, args.region_bed,
            args.blacklist_bed,
        )

    # WES: subtract blacklist from combined target+antitarget windows
    # (antitargets already subtract blacklist internally, but targets do not)
    # WGS: blacklist already handled inside generate_wgs_windows
    if args.blacklist_bed and wes_mode:
        n_before = len(windows)
        n_tgt_before = int(windows["is_target"].sum())
        windows = subtract_blacklist(windows, args.blacklist_bed)
        n_removed = n_before - len(windows)
        n_tgt_after = int(windows["is_target"].sum())
        n_tgt_rm = n_tgt_before - n_tgt_after
        n_anti_rm = n_removed - n_tgt_rm
        print(
            f"  {len(windows)} after blacklist subtraction "
            f"(removed {n_tgt_rm} target, {n_anti_rm} antitarget)"
        )

    # --- Step 2: Assign region_id, drop windows outside regions ---
    print("[2/6] Assigning region_id ...")
    windows["region_id"] = assign_region_id(windows, args.region_bed)
    n_before = len(windows)
    windows = windows[windows["region_id"].notna()].reset_index(drop=True)
    n_dropped = n_before - len(windows)
    assert n_dropped == 0, (
        f"{n_dropped}/{n_before} windows could not be assigned a region_id"
    )
    print(f"  {len(windows)} assigned")

    # --- Step 3: Compute GC ---
    print("[3/6] Computing GC content ...")
    windows = compute_gc(windows, args.reference)
    print(f"  GC range: [{windows['GC'].min():.4f}, {windows['GC'].max():.4f}]")
    out_cols = ["#CHR", "START", "END", "region_id", "GC"]

    # --- Step 4: Compute mappability (optional) ---
    if args.mappability_bed:
        print("[4/6] Computing mappability ...")
        windows["MAP"] = compute_mappability(
            windows, args.mappability_bed, args.genome_size,
        )
        out_cols.append("MAP")
    else:
        print("[4/6] Skipping mappability (no --mappability_bed)")

    # --- Step 5: Compute replication timing (optional) ---
    if args.repliseq:
        if args.reference_version not in ("hg38", "hg19"):
            print(f"[5/6] Skipping repli-seq (unsupported for {args.reference_version})")
        else:
            print("[5/6] Computing replication timing ...")
            windows["REPLI"] = compute_repliseq(windows, standard_chroms, args)
            out_cols.append("REPLI")
    else:
        print("[5/6] Skipping repli-seq (--repliseq not set)")

    if "is_target" in windows.columns:
        out_cols.append("is_target")

    # --- Summary ---
    print_summary(windows)

    # --- Step 6: Sort and write output ---
    windows = _sort_windows(windows)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    windows[out_cols].to_csv(
        args.out_file,
        sep="\t",
        header=True,
        index=False,
        compression="gzip" if args.out_file.endswith(".gz") else None,
    )
    elapsed = time.time() - t0
    print(f"[6/6] Wrote {len(windows)} rows to {args.out_file}")
    print(f"  Total elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
