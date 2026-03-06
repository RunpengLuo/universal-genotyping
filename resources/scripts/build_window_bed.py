#!/usr/bin/env python3
"""Build a per-window bias correction BED with region/blacklist/WES filtering.

Generates windows (WGS fixed-size or WES target+antitarget), filters by
accessible regions and blacklist, then computes per-window GC content with
optional mappability and replication timing columns.

WGS mode (default):
  1. Tile genome into fixed-size windows
  2. Keep windows overlapping --region_bed
  3. Subtract --blacklist_bed if provided

WES mode (--wes_targets_bed):
  1. Split vendor exon targets into ~target_avg_size pieces
  2. Generate antitarget bins in accessible non-target regions
  3. Combine target + antitarget, mark is_target column

Both modes then compute GC/MAP/REPLI and assign region_id.

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
    """Convert a DataFrame with #CHR/START/END to PyRanges with _idx."""
    win_pr = pr.PyRanges(
        chromosomes=df["#CHR"],
        starts=df["START"],
        ends=df["END"],
    )
    return win_pr.insert(pd.Series(np.arange(len(df)), name="_idx"))


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


# ---------------------------------------------------------------------------
# Reference version helpers
# ---------------------------------------------------------------------------


def get_standard_chroms(reference_version, genome_size_file):
    """Derive valid standard chromosomes from reference version and genome size file."""
    all_chroms = set(_read_chrom_sizes(genome_size_file)["chrom"])
    candidates = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
    if reference_version not in ("hg19", "hg38", "chm13v2"):
        # Fallback: also try unprefixed
        candidates |= {str(c) for c in list(range(1, 23)) + ["X", "Y"]}
    return candidates & all_chroms


# ---------------------------------------------------------------------------
# Window generation
# ---------------------------------------------------------------------------


def generate_windows(genome_size_file, window_size, standard_chroms):
    """Generate fixed-size windows from a chromosome sizes file."""
    sizes = _read_chrom_sizes(genome_size_file)
    sizes = sizes[sizes["chrom"].isin(standard_chroms)]

    rows = []
    for _, row in sizes.iterrows():
        for start in range(0, row["size"], window_size):
            rows.append((row["chrom"], start, min(start + window_size, row["size"])))

    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


def filter_windows_by_region(windows, region_bed_file):
    """Keep only windows overlapping accessible regions."""
    keep_idx = (
        _df_to_pyranges(windows)
        .overlap(pr.read_bed(region_bed_file))
        .df["_idx"]
        .unique()
    )
    return windows.iloc[keep_idx].sort_values(["#CHR", "START"]).reset_index(drop=True)


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
# WES window generation
# ---------------------------------------------------------------------------


def generate_wes_targets(wes_targets_bed, target_avg_size):
    """Subdivide large vendor exon targets into ~target_avg_size pieces."""
    df = pd.read_csv(
        wes_targets_bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["#CHR", "START", "END"],
        comment="#",
    )
    df = df[df["END"] > df["START"]].copy()
    rows = _subdivide_intervals(df, "#CHR", "START", "END", target_avg_size, min_size=1)
    return (
        pd.DataFrame(rows, columns=["#CHR", "START", "END"])
        .sort_values(["#CHR", "START"])
        .reset_index(drop=True)
    )


def generate_wes_antitargets(
    targets, region_bed_file, blacklist_bed_file, antitarget_avg_size, pad_size
):
    """Generate off-target bins in accessible regions not covered by exon targets."""
    access = pr.read_bed(region_bed_file)
    if blacklist_bed_file:
        access = access.subtract(pr.read_bed(blacklist_bed_file))

    tgt_pr = pr.PyRanges(
        chromosomes=targets["#CHR"],
        starts=targets["START"],
        ends=targets["END"],
    )
    padded = tgt_pr.copy()
    padded.Start = (padded.Start - pad_size).clip(lower=0)
    padded.End = padded.End + pad_size

    background = access.subtract(padded)
    background.Start = background.Start + pad_size
    background.End = background.End - pad_size
    bg_df = background.df
    bg_df = bg_df[bg_df["End"] > bg_df["Start"]].copy()

    min_bin_size = antitarget_avg_size // 16
    rows = _subdivide_intervals(
        bg_df, "Chromosome", "Start", "End", antitarget_avg_size, min_size=min_bin_size
    )
    return (
        pd.DataFrame(rows, columns=["#CHR", "START", "END"])
        .sort_values(["#CHR", "START"])
        .reset_index(drop=True)
    )


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
# GC content
# ---------------------------------------------------------------------------


def compute_gc(windows_df, reference):
    """Compute per-window GC fraction via pybedtools nucleotide_content."""
    bt = BedTool.from_dataframe(windows_df[["#CHR", "START", "END"]])
    nuc = bt.nucleotide_content(fi=reference).to_dataframe(disable_auto_names=True)
    windows_df = windows_df.copy()
    windows_df["GC"] = nuc["5_pct_gc"].values
    return windows_df


# ---------------------------------------------------------------------------
# Mappability
# ---------------------------------------------------------------------------


def compute_mappability(windows_df, mappability_file, genome_size_file):
    """Compute per-window mean mappability via pybedtools map."""
    bt = BedTool.from_dataframe(windows_df[["#CHR", "START", "END"]]).sort(
        g=genome_size_file
    )
    map_bt = BedTool(mappability_file).sort(g=genome_size_file)
    map_cov = bt.map(b=map_bt, c=4, o="mean", g=genome_size_file).to_dataframe(
        disable_auto_names=True
    )
    map_cov.columns = ["#CHR", "START", "END", "MAP"]
    return pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(0.0).clip(0.0, 1.0)


# ---------------------------------------------------------------------------
# Replication timing (Repli-seq)
# ---------------------------------------------------------------------------


def download(url, dest):
    """Download a file with wget."""
    subprocess.check_call(["wget", "-q", "-O", dest, url])


def bigwig_to_bedgraph(bigwig, bedgraph):
    """Convert bigWig to bedGraph using UCSC bigWigToBedGraph."""
    subprocess.check_call(["bigWigToBedGraph", bigwig, bedgraph])


def liftover(input_bed, chain, output_bed, unmapped):
    """Run UCSC liftOver."""
    subprocess.check_call(["liftOver", input_bed, chain, output_bed, unmapped])


def compute_repliseq(windows_df, standard_chroms, chain, bigwig_dir, work_dir):
    """Compute per-window mean replication timing from Repli-seq bigWigs.

    Downloads bigWig files if not present, converts to bedGraph, lifts from
    hg19 to hg38, and bins into the provided windows.
    """
    # --- Download bigWigs if needed ---
    os.makedirs(bigwig_dir, exist_ok=True)
    n_files = len(WAVE_SIGNAL_FILES)
    for i, fname in enumerate(WAVE_SIGNAL_FILES, 1):
        dest = os.path.join(bigwig_dir, fname)
        if not os.path.isfile(dest):
            print(f"\r  downloading bigWigs [{i}/{n_files}]", end="", flush=True)
            download(f"{UCSC_BASE}/{fname}", dest)
    print()

    bigwig_files = sorted(glob.glob(os.path.join(bigwig_dir, "*WaveSignal*.bigWig")))
    if not bigwig_files:
        sys.exit(f"Error: no WaveSignal bigWig files found in {bigwig_dir}")

    # --- Convert + liftOver each bigWig ---
    lifted_dir = os.path.join(work_dir, "lifted")
    os.makedirs(lifted_dir, exist_ok=True)
    n_bw = len(bigwig_files)

    lifted_files = []
    for i, bw in enumerate(bigwig_files, 1):
        name = os.path.basename(bw).replace(".bigWig", "")
        bg_file = os.path.join(lifted_dir, f"{name}.hg19.bedGraph")
        lifted_file = os.path.join(lifted_dir, f"{name}.hg38.bedGraph")
        unmapped_file = os.path.join(lifted_dir, f"{name}.unmapped")

        if not os.path.isfile(lifted_file):
            print(f"\r  liftOver [{i}/{n_bw}]", end="", flush=True)
            bigwig_to_bedgraph(bw, bg_file)
            liftover(bg_file, chain, lifted_file, unmapped_file)
            os.remove(bg_file)
        lifted_files.append(lifted_file)
    print()

    # --- Bin lifted signals into windows ---
    n_win = len(windows_df)
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int32)
    n_lf = len(lifted_files)

    for i, lf in enumerate(lifted_files, 1):
        print(f"\r  binning repli-seq [{i}/{n_lf}]", end="", flush=True)
        bg = pd.read_csv(
            lf,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "signal"],
            dtype={
                "chrom": str,
                "start": np.int64,
                "end": np.int64,
                "signal": np.float64,
            },
        )
        bg = bg[bg["chrom"].isin(standard_chroms)].reset_index(drop=True)

        for chrom, grp in bg.groupby("chrom", sort=False):
            win_mask = windows_df["#CHR"] == chrom
            win_idx = np.where(win_mask)[0]
            if len(win_idx) == 0:
                continue

            win_starts = windows_df["START"].to_numpy()[win_idx]
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
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
    parser.add_argument(
        "--reference_version",
        required=True,
        help="Reference genome version (hg19, hg38, chm13v2).",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference FASTA (must have .fai index).",
    )
    parser.add_argument(
        "--genome_size",
        required=True,
        help="Path to chromosome sizes file (tab-separated: chrom, size).",
    )
    parser.add_argument(
        "--region_bed",
        required=True,
        help="Accessible regions BED -- keep only windows overlapping these regions.",
    )
    parser.add_argument(
        "--out_file",
        required=True,
        help="Output gzipped TSV path (e.g. gc_map_repli.1kbp.hg38.bed.gz).",
    )
    parser.add_argument(
        "--blacklist_bed",
        default=None,
        help="Optional blacklist BED -- subtract these regions.",
    )
    parser.add_argument(
        "--wes_targets_bed",
        default=None,
        help=(
            "Vendor exon capture targets BED. If provided, enables WES mode "
            "(target+antitarget windows instead of fixed-size tiling)."
        ),
    )
    parser.add_argument(
        "--target_avg_size",
        type=int,
        default=200,
        help="Average target bin size for splitting large exon targets (WES only, default: 200).",
    )
    parser.add_argument(
        "--antitarget_avg_size",
        type=int,
        default=150000,
        help="Average off-target bin size (WES only, default: 150000).",
    )
    parser.add_argument(
        "--pad_size",
        type=int,
        default=500,
        help="Padding around targets for antitarget generation (WES only, default: 500).",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=1000,
        help="Window size in bp (WGS only, default: 1000).",
    )
    parser.add_argument(
        "--mappability_bed",
        default=None,
        help="Optional BED-format mappability track (4th column = score).",
    )
    parser.add_argument(
        "--repliseq",
        action="store_true",
        help="Enable replication timing from ENCODE Repli-seq bigWigs.",
    )
    parser.add_argument(
        "--chain",
        default=None,
        help=(
            "hg19-to-hg38 liftOver chain file. "
            "If not provided and --repliseq is set, downloads hg19ToHg38.over.chain.gz."
        ),
    )
    parser.add_argument(
        "--bigwig_dir",
        default=None,
        help=(
            "Directory containing pre-downloaded WaveSignal bigWig files. "
            "If not provided, bigWigs are downloaded to --work_dir/bigwigs/."
        ),
    )
    parser.add_argument(
        "--work_dir",
        default=None,
        help="Working directory for intermediate files (default: temp dir).",
    )
    args = parser.parse_args()

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
    if not wes_mode:
        print(f"  Window:      {args.window} bp")
    else:
        print(f"  WES targets: {args.wes_targets_bed}")
        print(f"  Target avg:  {args.target_avg_size} bp")
        print(f"  Antitarget:  {args.antitarget_avg_size} bp")
        print(f"  Pad size:    {args.pad_size} bp")
    print(f"  Mappability: {args.mappability_bed or 'none'}")
    print(f"  Repli-seq:   {'yes' if args.repliseq else 'no'}")
    print(f"  Output:      {args.out_file}")
    print()

    standard_chroms = get_standard_chroms(args.reference_version, args.genome_size)
    print(f"  {len(standard_chroms)} standard chromosomes")

    # --- Step 1: Generate windows ---
    if wes_mode:
        print("[1/6] WES: generating target + antitarget windows ...")
        targets = generate_wes_targets(args.wes_targets_bed, args.target_avg_size)
        antitargets = generate_wes_antitargets(
            targets,
            args.region_bed,
            args.blacklist_bed,
            args.antitarget_avg_size,
            args.pad_size,
        )
        targets["is_target"] = 1
        antitargets["is_target"] = 0
        windows = pd.concat([targets, antitargets], ignore_index=True)
        windows = windows.sort_values(["#CHR", "START"]).reset_index(drop=True)
        windows = windows[windows["#CHR"].isin(standard_chroms)].reset_index(drop=True)
        n_tgt = int(windows["is_target"].sum())
        print(
            f"  {len(windows)} windows ({n_tgt} target, {len(windows) - n_tgt} antitarget)"
        )
    else:
        print("[1/6] WGS: generating fixed windows ...")
        windows = generate_windows(args.genome_size, args.window, standard_chroms)
        print(
            f"  {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes"
        )

        windows = filter_windows_by_region(windows, args.region_bed)
        print(f"  {len(windows)} after region filter")

        if args.blacklist_bed:
            windows = subtract_blacklist(windows, args.blacklist_bed)
            print(f"  {len(windows)} after blacklist subtraction")

    # --- Step 2: Assign region_id ---
    print("[2/6] Assigning region_id ...")
    windows["region_id"] = assign_region_id(windows, args.region_bed)
    n_assigned = int(windows["region_id"].notna().sum())
    print(f"  {n_assigned}/{len(windows)} assigned")

    # --- Step 3: Compute GC ---
    print("[3/6] Computing GC content ...")
    windows = compute_gc(windows, args.reference)
    print(f"  GC range: [{windows['GC'].min():.4f}, {windows['GC'].max():.4f}]")
    out_cols = ["#CHR", "START", "END", "region_id", "GC"]

    # --- Step 4: Compute mappability (optional) ---
    if args.mappability_bed:
        print("[4/6] Computing mappability ...")
        windows["MAP"] = compute_mappability(
            windows, args.mappability_bed, args.genome_size
        )
        out_cols.append("MAP")
    else:
        print("[4/6] Skipping mappability (no --mappability_bed)")

    # --- Step 5: Compute replication timing (optional) ---
    if args.repliseq:
        if args.reference_version not in ("hg38", "hg19"):
            print(
                f"[5/6] Skipping repli-seq (unsupported for {args.reference_version})"
            )
        else:
            print("[5/6] Computing replication timing ...")
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
                    download(CHAIN_URL, chain)
            if not os.path.isfile(chain):
                sys.exit(f"Error: chain file not found: {chain}")

            bigwig_dir = args.bigwig_dir or os.path.join(work_dir, "bigwigs")
            windows["REPLI"] = compute_repliseq(
                windows, standard_chroms, chain, bigwig_dir, work_dir
            )
            out_cols.append("REPLI")

            if cleanup:
                import shutil

                shutil.rmtree(work_dir)
    else:
        print("[5/6] Skipping repli-seq (--repliseq not set)")

    if "is_target" in windows.columns:
        out_cols.append("is_target")

    # --- Summary stats ---
    sizes = windows["END"] - windows["START"]
    n_regions = windows["region_id"].nunique()
    n_chroms = windows["#CHR"].nunique()
    print()
    print("  Summary:")
    print(f"    {len(windows)} windows across {n_chroms} chromosomes")
    print(f"    {n_regions} distinct regions")
    print(
        f"    window size: min={sizes.min()}, median={int(sizes.median())}, "
        f"max={sizes.max()}, mean={sizes.mean():.0f}"
    )
    total_bp = int(sizes.sum())
    print(f"    total coverage: {total_bp:,} bp ({total_bp / 1e9:.2f} Gb)")
    if "is_target" in windows.columns:
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

    # Per-chromosome stats
    has_target = "is_target" in windows.columns
    header = f"    {'chrom':<8} {'windows':>8} {'regions':>8} {'coverage_Mb':>12} {'median_size':>12}"
    if has_target:
        header += f" {'target':>8} {'antitarget':>10}"
    print(header)
    # Sort chromosomes naturally (chr1, chr2, ..., chr22, chrX, chrY)
    chrom_order = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    for chrom in chrom_order:
        if chrom not in windows["#CHR"].values:
            continue
        mask = windows["#CHR"] == chrom
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

    # --- Step 6: Write output ---
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
