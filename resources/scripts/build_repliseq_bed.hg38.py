#!/usr/bin/env python3
"""LiftOver UCSC ENCODE/UW Repli-seq WaveSignal bigWigs from hg19 to hg38
and build a per-window consensus replication timing BED.

Downloads wavelet-smoothed Repli-seq signal bigWig files (hg19) from:
  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/

For each cell line the script converts to bedGraph, lifts over to hg38,
then averages across all cell lines into fixed-size genomic windows.

Output is a gzipped TSV: #CHR  START  END  REPLI

Requirements:
  - bigWigToBedGraph  (UCSC tools)
  - liftOver          (UCSC tools)
  - hg19ToHg38.over.chain.gz  (from UCSC goldenPath)
  - pandas, numpy

Cell lines with WaveSignal bigWig (16 total):
  Bg02es, Bj, Gm06990, Gm12801, Gm12812, Gm12813, Gm12878,
  Helas3, Hepg2, Huvec, Imr90, K562, Mcf7, Nhek, Sknsh
  (Bj has two replicates; all others have one)

Reference:
  Hansen et al. (2010) Genome-wide replication profiles from ENCODE.
  https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeUwRepliSeq
"""

import argparse
import glob
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


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


def download(url, dest):
    """Download a file with wget."""
    subprocess.check_call(["wget", "-q", "-O", dest, url])


def bigwig_to_bedgraph(bigwig, bedgraph):
    """Convert bigWig to bedGraph using UCSC bigWigToBedGraph."""
    subprocess.check_call(["bigWigToBedGraph", bigwig, bedgraph])


def liftover(input_bed, chain, output_bed, unmapped):
    """Run UCSC liftOver."""
    subprocess.check_call(["liftOver", input_bed, chain, output_bed, unmapped])


def generate_windows(genome_size_file, window_size):
    """Generate fixed-size windows from a chromosome sizes file."""
    sizes = pd.read_csv(
        genome_size_file, sep="\t", header=None, names=["chrom", "size"],
        dtype={"chrom": str, "size": np.int64},
    )
    sizes = sizes[sizes["chrom"].isin(STANDARD_CHROMS)]
    rows = []
    for _, row in sizes.iterrows():
        for start in range(0, row["size"], window_size):
            rows.append((row["chrom"], start, min(start + window_size, row["size"])))
    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a consensus replication timing BED (hg38) from UCSC "
            "ENCODE/UW Repli-seq WaveSignal bigWigs (hg19)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python build_repliseq_bed.hg38.py \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --chain /path/to/hg19ToHg38.over.chain.gz \\\n"
            "    --out_file repliseq.1kbp.hg38.bed.gz\n\n"
            "The script will download bigWig files to --work_dir on first run.\n"
            "Re-running with the same --work_dir skips already-downloaded files."
        ),
    )
    parser.add_argument(
        "--genome_size", required=True,
        help="hg38 chromosome sizes file (tab-separated: chrom, size).",
    )
    parser.add_argument(
        "--chain", default=None,
        help=(
            "hg19-to-hg38 liftOver chain file. "
            "If not provided, downloads hg19ToHg38.over.chain.gz to --work_dir."
        ),
    )
    parser.add_argument(
        "--out_file", required=True,
        help="Output gzipped TSV (e.g. repliseq.1kbp.hg38.bed.gz).",
    )
    parser.add_argument(
        "--window", type=int, default=1000,
        help="Window size in bp (default: 1000).",
    )
    parser.add_argument(
        "--work_dir", default=None,
        help="Working directory for intermediate files (default: temp dir).",
    )
    parser.add_argument(
        "--bigwig_dir", default=None,
        help=(
            "Directory containing pre-downloaded WaveSignal bigWig files. "
            "If not provided, bigWigs are downloaded to --work_dir/bigwigs/."
        ),
    )
    args = parser.parse_args()

    if not os.path.isfile(args.genome_size):
        sys.exit(f"Error: genome_size not found: {args.genome_size}")

    # Set up working directory
    if args.work_dir:
        work_dir = args.work_dir
        os.makedirs(work_dir, exist_ok=True)
        cleanup = False
    else:
        work_dir = tempfile.mkdtemp(prefix="repliseq_")
        cleanup = True

    print("=== build_repliseq_bed.hg38.py ===")
    print(f"  Genome size: {args.genome_size}")
    print(f"  Window:      {args.window} bp")
    print(f"  Work dir:    {work_dir}")
    print(f"  Output:      {args.out_file}")
    print()

    # --- Step 1: Obtain chain file ---
    chain = args.chain
    if chain is None:
        chain = os.path.join(work_dir, "hg19ToHg38.over.chain.gz")
        if not os.path.isfile(chain):
            print("[Step 1] Downloading chain file ...")
            download(CHAIN_URL, chain)
        else:
            print("[Step 1] Chain file already present.")
    else:
        print(f"[Step 1] Using provided chain: {chain}")
    if not os.path.isfile(chain):
        sys.exit(f"Error: chain file not found: {chain}")

    # --- Step 2: Download / locate bigWig files ---
    bigwig_dir = args.bigwig_dir
    if bigwig_dir is None:
        bigwig_dir = os.path.join(work_dir, "bigwigs")
        os.makedirs(bigwig_dir, exist_ok=True)
        print(f"[Step 2] Downloading bigWig files to {bigwig_dir} ...")
        for fname in WAVE_SIGNAL_FILES:
            dest = os.path.join(bigwig_dir, fname)
            if os.path.isfile(dest):
                print(f"  skip (exists): {fname}")
                continue
            print(f"  downloading: {fname}")
            download(f"{UCSC_BASE}/{fname}", dest)
    else:
        print(f"[Step 2] Using bigWigs from {bigwig_dir}")

    bigwig_files = sorted(glob.glob(os.path.join(bigwig_dir, "*WaveSignal*.bigWig")))
    if not bigwig_files:
        sys.exit(f"Error: no WaveSignal bigWig files found in {bigwig_dir}")
    print(f"  Found {len(bigwig_files)} bigWig files")

    # --- Step 3: Convert + liftOver each bigWig ---
    lifted_dir = os.path.join(work_dir, "lifted")
    os.makedirs(lifted_dir, exist_ok=True)
    print(f"[Step 3] Converting bigWig → bedGraph → liftOver hg38 ...")

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
        # clean up hg19 bedGraph to save space
        os.remove(bg_file)
        lifted_files.append(lifted_file)

    # --- Step 4: Generate windows and compute per-window mean ---
    print(f"[Step 4] Generating {args.window} bp windows ...")
    windows = generate_windows(args.genome_size, args.window)
    n_win = len(windows)
    print(f"  {n_win} windows across {windows['#CHR'].nunique()} chromosomes")

    # Accumulate signal sums and counts per window
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int32)

    print(f"[Step 5] Binning lifted signals into windows ...")
    for lf in lifted_files:
        name = os.path.basename(lf)
        bg = pd.read_csv(
            lf, sep="\t", header=None,
            names=["chrom", "start", "end", "signal"],
            dtype={"chrom": str, "start": np.int64, "end": np.int64, "signal": np.float64},
        )
        bg = bg[bg["chrom"].isin(STANDARD_CHROMS)].reset_index(drop=True)

        for chrom, grp in bg.groupby("chrom", sort=False):
            win_mask = windows["#CHR"] == chrom
            win_idx = np.where(win_mask)[0]
            if len(win_idx) == 0:
                continue

            win_starts = windows["START"].to_numpy()[win_idx]
            # assign each bedGraph interval to its overlapping window
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

    windows["REPLI"] = np.round(mean_signal, 6)

    # --- Step 6: Write output ---
    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    windows.to_csv(
        args.out_file, sep="\t", header=True, index=False,
        compression="gzip" if args.out_file.endswith(".gz") else None,
    )
    print(f"[Step 6] Wrote {n_win} rows to {args.out_file}")

    if cleanup:
        import shutil
        shutil.rmtree(work_dir)
        print(f"  cleaned up temp dir: {work_dir}")


if __name__ == "__main__":
    main()
