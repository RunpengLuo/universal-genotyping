#!/usr/bin/env python3
"""Build a fixed-window bias correction BED for hg38.

Computes per-window GC content from a reference FASTA, with optional
mappability and replication timing columns.

GC and mappability are computed via pybedtools. Replication timing is
derived from UCSC ENCODE/UW Repli-seq WaveSignal bigWigs (hg19),
lifted over to hg38 via UCSC liftOver.

Output is a gzipped TSV: #CHR  START  END  GC  [MAP]  [REPLI]

Dependencies:
  - pybedtools, pandas, numpy
  - bigWigToBedGraph, liftOver (UCSC tools) — only if --repliseq
"""

import argparse
import glob
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from pybedtools import BedTool


STANDARD_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}

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
    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/"
    "hg19ToHg38.over.chain.gz"
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def generate_windows(genome_size_file, window_size):
    """Generate fixed-size windows from a chromosome sizes file.

    Parameters
    ----------
    genome_size_file : str
        Tab-separated file with columns (chrom, size).
    window_size : int
        Window size in base pairs.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``#CHR``, ``START``, ``END``.
    """
    sizes = pd.read_csv(
        genome_size_file, sep="\t", header=None, names=["chrom", "size"],
        dtype={"chrom": str, "size": np.int64},
    )
    sizes = sizes[sizes["chrom"].isin(STANDARD_CHROMS)]

    rows = []
    for _, row in sizes.iterrows():
        chrom = row["chrom"]
        chrom_len = row["size"]
        for start in range(0, chrom_len, window_size):
            end = min(start + window_size, chrom_len)
            rows.append((chrom, start, end))

    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


# ---------------------------------------------------------------------------
# GC content
# ---------------------------------------------------------------------------

def compute_gc(windows_df, reference):
    """Compute per-window GC fraction via pybedtools nucleotide_content.

    Parameters
    ----------
    windows_df : pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` columns.
    reference : str
        Path to the reference FASTA (must be indexed).

    Returns
    -------
    pd.DataFrame
        Input DataFrame with an added ``GC`` column.
    """
    bt = BedTool.from_dataframe(windows_df[["#CHR", "START", "END"]])
    nuc = bt.nucleotide_content(fi=reference).to_dataframe(disable_auto_names=True)
    nuc = nuc.rename(columns={
        "#1_usercol": "#CHR",
        "2_usercol": "START",
        "3_usercol": "END",
        "5_pct_gc": "GC",
    })
    result = windows_df.merge(nuc[["#CHR", "START", "END", "GC"]])
    return result


# ---------------------------------------------------------------------------
# Mappability
# ---------------------------------------------------------------------------

def compute_mappability(windows_df, mappability_file, genome_size_file):
    """Compute per-window mean mappability via pybedtools map.

    Parameters
    ----------
    windows_df : pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` columns.
    mappability_file : str
        Path to a BED-format mappability track (4th column = score).
    genome_size_file : str
        Path to chromosome sizes file.

    Returns
    -------
    pd.Series
        Per-window mean mappability values.
    """
    bt = BedTool.from_dataframe(windows_df[["#CHR", "START", "END"]]).sort(g=genome_size_file)
    map_bt = BedTool(mappability_file).sort(g=genome_size_file)
    map_cov = bt.map(b=map_bt, c=4, o="mean", g=genome_size_file).to_dataframe(
        disable_auto_names=True
    )
    map_cov.columns = ["#CHR", "START", "END", "MAP"]
    map_vals = pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    return map_vals


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


def compute_repliseq(windows_df, chain, bigwig_dir, work_dir):
    """Compute per-window mean replication timing from Repli-seq bigWigs.

    Downloads bigWig files if not present, converts to bedGraph, lifts from
    hg19 to hg38, and bins into the provided windows.

    Parameters
    ----------
    windows_df : pd.DataFrame
        DataFrame with ``#CHR``, ``START``, ``END`` columns.
    chain : str
        Path to hg19→hg38 liftOver chain file.
    bigwig_dir : str
        Directory containing (or to download) WaveSignal bigWig files.
    work_dir : str
        Directory for intermediate files (bedGraph, lifted).

    Returns
    -------
    np.ndarray
        Per-window mean replication timing signal (NaN where no data).
    """
    # --- Download bigWigs if needed ---
    os.makedirs(bigwig_dir, exist_ok=True)
    for fname in WAVE_SIGNAL_FILES:
        dest = os.path.join(bigwig_dir, fname)
        if os.path.isfile(dest):
            print(f"  skip (exists): {fname}")
            continue
        print(f"  downloading: {fname}")
        download(f"{UCSC_BASE}/{fname}", dest)

    bigwig_files = sorted(glob.glob(os.path.join(bigwig_dir, "*WaveSignal*.bigWig")))
    if not bigwig_files:
        sys.exit(f"Error: no WaveSignal bigWig files found in {bigwig_dir}")
    print(f"  Found {len(bigwig_files)} bigWig files")

    # --- Convert + liftOver each bigWig ---
    lifted_dir = os.path.join(work_dir, "lifted")
    os.makedirs(lifted_dir, exist_ok=True)
    print("  Converting bigWig -> bedGraph -> liftOver hg38 ...")

    lifted_files = []
    for bw in bigwig_files:
        name = os.path.basename(bw).replace(".bigWig", "")
        bg_file = os.path.join(lifted_dir, f"{name}.hg19.bedGraph")
        lifted_file = os.path.join(lifted_dir, f"{name}.hg38.bedGraph")
        unmapped_file = os.path.join(lifted_dir, f"{name}.unmapped")

        if os.path.isfile(lifted_file):
            print(f"  skip (exists): {name}")
            lifted_files.append(lifted_file)
            continue

        print(f"  processing: {name}")
        bigwig_to_bedgraph(bw, bg_file)
        liftover(bg_file, chain, lifted_file, unmapped_file)
        os.remove(bg_file)
        lifted_files.append(lifted_file)

    # --- Bin lifted signals into windows ---
    n_win = len(windows_df)
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int32)

    print("  Binning lifted signals into windows ...")
    for lf in lifted_files:
        name = os.path.basename(lf)
        bg = pd.read_csv(
            lf, sep="\t", header=None,
            names=["chrom", "start", "end", "signal"],
            dtype={"chrom": str, "start": np.int64, "end": np.int64, "signal": np.float64},
        )
        bg = bg[bg["chrom"].isin(STANDARD_CHROMS)].reset_index(drop=True)

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

        print(f"  binned: {name}")

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
            "Build a fixed-window bias correction BED for hg38. "
            "Computes per-window GC content, with optional mappability "
            "and replication timing columns."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # GC-only (1 kb windows)\n"
            "  python build_bias_bed.hg38.py \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --out_file gc.1kbp.hg38.bed.gz\n\n"
            "  # GC + mappability\n"
            "  python build_bias_bed.hg38.py \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --mappability_bed /path/to/mappability.bed \\\n"
            "    --out_file gc_map.1kbp.hg38.bed.gz\n\n"
            "  # GC + mappability + replication timing\n"
            "  python build_bias_bed.hg38.py \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --mappability_bed /path/to/mappability.bed \\\n"
            "    --repliseq \\\n"
            "    --out_file gc_map_repli.1kbp.hg38.bed.gz"
        ),
    )
    parser.add_argument(
        "--reference", required=True,
        help="Path to reference FASTA (must have .fai index).",
    )
    parser.add_argument(
        "--genome_size", required=True,
        help="Path to chromosome sizes file (tab-separated: chrom, size).",
    )
    parser.add_argument(
        "--out_file", required=True,
        help="Output gzipped TSV path (e.g. gc_map_repli.1kbp.hg38.bed.gz).",
    )
    parser.add_argument(
        "--window", type=int, default=1000,
        help="Window size in bp (default: 1000).",
    )
    parser.add_argument(
        "--mappability_bed", default=None,
        help="Optional BED-format mappability track (4th column = score).",
    )
    parser.add_argument(
        "--repliseq", action="store_true",
        help="Enable replication timing from ENCODE Repli-seq bigWigs.",
    )
    parser.add_argument(
        "--chain", default=None,
        help=(
            "hg19-to-hg38 liftOver chain file. "
            "If not provided and --repliseq is set, downloads hg19ToHg38.over.chain.gz."
        ),
    )
    parser.add_argument(
        "--bigwig_dir", default=None,
        help=(
            "Directory containing pre-downloaded WaveSignal bigWig files. "
            "If not provided, bigWigs are downloaded to --work_dir/bigwigs/."
        ),
    )
    parser.add_argument(
        "--work_dir", default=None,
        help="Working directory for intermediate files (default: temp dir).",
    )
    args = parser.parse_args()

    # --- Validate inputs ---
    if not os.path.isfile(args.reference):
        sys.exit(f"Error: reference not found: {args.reference}")
    if not os.path.isfile(args.genome_size):
        sys.exit(f"Error: genome_size not found: {args.genome_size}")
    if args.mappability_bed and not os.path.isfile(args.mappability_bed):
        sys.exit(f"Error: mappability_bed not found: {args.mappability_bed}")

    print("=== build_bias_bed.hg38.py ===")
    print(f"  Reference:   {args.reference}")
    print(f"  Genome size: {args.genome_size}")
    print(f"  Window:      {args.window} bp")
    print(f"  Mappability: {args.mappability_bed or 'none (MAP column omitted)'}")
    print(f"  Repli-seq:   {'yes' if args.repliseq else 'no (REPLI column omitted)'}")
    print(f"  Output:      {args.out_file}")
    print()

    # --- Step 1: Generate windows ---
    print("[Step 1] Generating fixed windows ...")
    windows = generate_windows(args.genome_size, args.window)
    print(f"  {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes")

    # --- Step 2: Compute GC ---
    print("[Step 2] Computing GC content ...")
    result = compute_gc(windows, args.reference)
    print(f"  GC range: [{result['GC'].min():.4f}, {result['GC'].max():.4f}]")
    out_cols = ["#CHR", "START", "END", "GC"]

    # --- Step 3: Compute mappability (optional) ---
    if args.mappability_bed:
        print("[Step 3] Computing mappability ...")
        result["MAP"] = compute_mappability(result, args.mappability_bed, args.genome_size)
        out_cols.append("MAP")
    else:
        print("[Step 3] No mappability track; skipping MAP column")

    # --- Step 4: Compute replication timing (optional) ---
    if args.repliseq:
        # Set up working directory
        if args.work_dir:
            work_dir = args.work_dir
            os.makedirs(work_dir, exist_ok=True)
            cleanup = False
        else:
            work_dir = tempfile.mkdtemp(prefix="repliseq_")
            cleanup = True

        # Obtain chain file
        chain = args.chain
        if chain is None:
            chain = os.path.join(work_dir, "hg19ToHg38.over.chain.gz")
            if not os.path.isfile(chain):
                print("[Step 4a] Downloading chain file ...")
                download(CHAIN_URL, chain)
            else:
                print("[Step 4a] Chain file already present.")
        else:
            print(f"[Step 4a] Using provided chain: {chain}")
        if not os.path.isfile(chain):
            sys.exit(f"Error: chain file not found: {chain}")

        bigwig_dir = args.bigwig_dir or os.path.join(work_dir, "bigwigs")

        print("[Step 4b] Computing replication timing ...")
        result["REPLI"] = compute_repliseq(result, chain, bigwig_dir, work_dir)
        out_cols.append("REPLI")

        if cleanup:
            import shutil
            shutil.rmtree(work_dir)
            print(f"  cleaned up temp dir: {work_dir}")
    else:
        print("[Step 4] Repli-seq disabled; skipping REPLI column")

    # --- Step 5: Write output ---
    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    result[out_cols].to_csv(
        args.out_file, sep="\t", header=True, index=False,
        compression="gzip" if args.out_file.endswith(".gz") else None,
    )
    print(f"[Step 5] Wrote {len(result)} rows to {args.out_file}")


if __name__ == "__main__":
    main()
