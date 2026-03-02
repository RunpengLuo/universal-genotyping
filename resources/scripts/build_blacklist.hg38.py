#!/usr/bin/env python3
"""Build a merged repeat/segdup/low-mappability blacklist BED for hg38.

Combines the functionality of the former build_mappability.sh and
build_blacklist.hg38.py into a single script.  Given a reference FASTA,
it runs hmmcopy_utils (generateMap.pl + mapCounter) to produce per-base
mappability, bins and filters low-mappability regions, downloads UCSC
repeat/segdup tables, merges everything with pyranges, and writes a
gzipped blacklist BED.

Dependencies on PATH:
  - generateMap.pl, mapCounter, bowtie  (hmmcopy_utils + bowtie)
  - bowtie-build  (only if --build_index is used)
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request

import numpy as np
import pandas as pd
import pyranges as pr

UCSC_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"
TABLES = ["genomicSuperDups.txt.gz", "simpleRepeat.txt.gz"]
STANDARD_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def _reporthook(block_num, block_size, total_size):
    """Print download progress."""
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(downloaded / total_size * 100, 100)
        mb = downloaded / 1e6
        total_mb = total_size / 1e6
        sys.stdout.write(f"\r  {mb:.1f}/{total_mb:.1f} MB ({pct:.0f}%)")
    else:
        sys.stdout.write(f"\r  {downloaded / 1e6:.1f} MB")
    sys.stdout.flush()


def download_file(url, dest_dir, filename=None):
    """Download a file to *dest_dir* and return the local path."""
    if filename is None:
        filename = os.path.basename(url)
    local_path = os.path.join(dest_dir, filename)
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, local_path, reporthook=_reporthook)
    print()  # newline after progress
    return local_path


def read_bed_columns(path):
    """Read chrom/start/end (columns 1-3) from a gzipped UCSC table file."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[1, 2, 3],
        names=["Chromosome", "Start", "End"],
        dtype={"Chromosome": str, "Start": int, "End": int},
    )
    df = df[df["Chromosome"].isin(STANDARD_CHROMS)]
    return df


def check_dependencies(build_index):
    """Verify required executables are on PATH."""
    required = ["generateMap.pl", "mapCounter", "bowtie"]
    if build_index:
        required.append("bowtie-build")
    missing = [cmd for cmd in required if shutil.which(cmd) is None]
    if missing:
        sys.exit(
            f"Error: the following commands were not found on PATH: {', '.join(missing)}\n"
            "Install hmmcopy_utils (https://github.com/shahcompbio/hmmcopy_utils) "
            "and bowtie (https://bowtie-bio.sourceforge.net)."
        )


def run_bowtie_build(reference, index_base, threads):
    """Build a bowtie index from the reference FASTA."""
    print(f"[Step 0] Building bowtie index at {index_base} ...")
    subprocess.run(
        ["bowtie-build", reference, index_base, "--threads", str(threads)],
        check=True,
    )
    print("  Done.")


def run_generate_map(reference, index_base, map_bw, read_length, threads):
    """Run generateMap.pl to produce a per-base mappability BigWig."""
    print(f"[Step 1] Running generateMap.pl (read_length={read_length}) ...")
    subprocess.run(
        [
            "generateMap.pl",
            "-b", index_base,
            reference,
            "-l", str(read_length),
            "-o", map_bw,
            "-t", str(threads),
        ],
        check=True,
    )
    print(f"  Output: {map_bw}")


def run_map_counter(bw_path, window, chroms):
    """Run hmmcopy_utils mapCounter on a per-base mappability BigWig and return
    a DataFrame of bins with their average mappability scores.

    mapCounter outputs fixedStep WIG:
        fixedStep chrom=chr1 start=1 step=1000 span=1000
        0.235191
        0.905116
        ...

    Returns a DataFrame with columns [Chromosome, Start, End, Score].
    """
    chrom_str = ",".join(sorted(chroms, key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c[3:])))
    cmd = ["mapCounter", "-w", str(window), "-c", chrom_str, bw_path]
    print(f"[Step 2] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    rows = []
    chrom = None
    start = None
    step = None
    span = None

    for line in result.stdout.splitlines():
        if line.startswith("fixedStep"):
            m = re.match(
                r"fixedStep\s+chrom=(\S+)\s+start=(\d+)\s+step=(\d+)\s+span=(\d+)",
                line,
            )
            if m is None:
                continue
            chrom = m.group(1)
            start = int(m.group(2)) - 1  # WIG is 1-based, convert to 0-based
            step = int(m.group(3))
            span = int(m.group(4))
        elif chrom is not None and line.strip():
            score = float(line.strip())
            rows.append((chrom, start, start + span, score))
            start += step

    if not rows:
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "Score"])

    arr = np.array(rows, dtype=object)
    return pd.DataFrame({
        "Chromosome": arr[:, 0],
        "Start": arr[:, 1].astype(np.int64),
        "End": arr[:, 2].astype(np.int64),
        "Score": arr[:, 3].astype(np.float64),
    })


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a merged repeat/segdup/low-mappability blacklist BED for hg38. "
            "Runs hmmcopy_utils (generateMap.pl + mapCounter) to compute mappability "
            "from a reference FASTA, downloads UCSC repeat tables, and merges "
            "everything into a single gzipped blacklist BED."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Dependencies:\n"
            "  generateMap.pl, mapCounter, bowtie  (hmmcopy_utils + bowtie)\n"
            "  bowtie-build  (only with --build_index)\n\n"
            "Example:\n"
            "  python build_blacklist.hg38.py \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --out_file /path/to/blacklist.hg38.bed.gz \\\n"
            "    --build_index --threads 8"
        ),
    )
    parser.add_argument(
        "--reference", required=True,
        help="Path to reference FASTA (must be indexed: .fai).",
    )
    parser.add_argument(
        "--out_file", required=True,
        help="Output gzipped BED file path (e.g. blacklist.hg38.bed.gz).",
    )
    parser.add_argument(
        "--window", type=int, default=1000,
        help="Bin size (bp) for mapCounter (default: 1000).",
    )
    parser.add_argument(
        "--read_length", type=int, default=150,
        help="Read length for generateMap.pl (default: 150).",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for bowtie alignment (default: 4).",
    )
    parser.add_argument(
        "--min_map_score", type=float, default=0.9,
        help="Exclude bins with average mappability below this threshold (default: 0.9).",
    )
    parser.add_argument(
        "--build_index", action="store_true",
        help="Build a bowtie index before running generateMap.pl.",
    )
    parser.add_argument(
        "--work_dir", default=None,
        help="Directory for intermediate files (default: auto temp dir, cleaned up on exit).",
    )
    args = parser.parse_args()

    # Validate reference
    if not os.path.isfile(args.reference):
        sys.exit(f"Error: reference file not found: {args.reference}")

    reference = os.path.abspath(args.reference)

    # Check dependencies
    check_dependencies(args.build_index)

    # Set up work directory
    user_work_dir = args.work_dir is not None
    if user_work_dir:
        work_dir = os.path.abspath(args.work_dir)
        os.makedirs(work_dir, exist_ok=True)
    else:
        work_dir = tempfile.mkdtemp(prefix="blacklist_")

    print("=== build_blacklist.hg38.py ===")
    print(f"  Reference:     {reference}")
    print(f"  Output:        {args.out_file}")
    print(f"  Window:        {args.window}")
    print(f"  Read length:   {args.read_length}")
    print(f"  Threads:       {args.threads}")
    print(f"  Min map score: {args.min_map_score}")
    print(f"  Build index:   {args.build_index}")
    print(f"  Work dir:      {work_dir}")
    print()

    beds = []
    try:
        index_base = os.path.join(work_dir, "ref_index")
        map_bw = os.path.join(work_dir, "mappability.bw")

        # Step 0 (optional): build bowtie index
        if args.build_index:
            run_bowtie_build(reference, index_base, args.threads)
        else:
            if not (
                os.path.isfile(f"{index_base}.1.ebwt")
                or os.path.isfile(f"{index_base}.1.ebwtl")
            ):
                print(
                    f"Warning: no bowtie index found at {index_base}.*\n"
                    "  Use --build_index to create one, or ensure generateMap.pl\n"
                    "  can locate the index via its own defaults."
                )

        # Step 1: generateMap.pl → per-base BigWig
        run_generate_map(reference, index_base, map_bw, args.read_length, args.threads)

        # Step 2: mapCounter → binned DataFrame
        binned = run_map_counter(map_bw, args.window, STANDARD_CHROMS)

        # Step 3: filter low-mappability bins
        low_map = binned[binned["Score"] < args.min_map_score][["Chromosome", "Start", "End"]].copy()
        total_bins = len(binned)
        n_low = len(low_map)
        pct = n_low / total_bins * 100 if total_bins > 0 else 0
        print(
            f"[Step 3] Low-mappability bins (score < {args.min_map_score}): "
            f"{n_low}/{total_bins} ({pct:.1f}%)"
        )
        beds.append(low_map)

        # Step 4: UCSC repeat/segdup tables
        print("[Step 4] Downloading UCSC repeat tables ...")
        for table in TABLES:
            local = download_file(f"{UCSC_BASE}/{table}", work_dir, table)
            bed = read_bed_columns(local)
            print(f"  {table}: {len(bed)} intervals")
            beds.append(bed)

    finally:
        if not user_work_dir:
            shutil.rmtree(work_dir, ignore_errors=True)

    # Step 5: merge and write
    combined = pd.concat(beds, ignore_index=True)
    merged = pr.PyRanges(combined).merge().sort()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    merged.as_df()[["Chromosome", "Start", "End"]].to_csv(
        args.out_file, sep="\t", header=False, index=False,
        compression="gzip" if args.out_file.endswith(".gz") else None,
    )
    print(f"[Step 5] Wrote {len(merged)} intervals to {args.out_file}")


if __name__ == "__main__":
    main()
