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

NCBI_REFERENCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)
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
    """Download a file to *dest_dir* and return the local path.

    Skips the download if the file already exists.
    """
    if filename is None:
        filename = os.path.basename(url)
    local_path = os.path.join(dest_dir, filename)
    if os.path.isfile(local_path):
        print(f"  [skip] {local_path} already exists")
        return local_path
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


def download_reference(dest_dir):
    """Download the NCBI GRCh38 analysis set FASTA, decompress, and index.

    Skips steps whose outputs already exist.
    """
    gz_path = os.path.join(dest_dir, os.path.basename(NCBI_REFERENCE_URL))
    fa_path = gz_path.removesuffix(".gz")
    if os.path.isfile(fa_path):
        print(f"  [skip] Reference already exists: {fa_path}")
    else:
        download_file(NCBI_REFERENCE_URL, dest_dir)
        print("  Decompressing ...")
        subprocess.run(["gunzip", gz_path], check=True)
    if os.path.isfile(fa_path + ".fai"):
        print(f"  [skip] Index already exists: {fa_path}.fai")
    else:
        print("  Indexing with samtools faidx ...")
        subprocess.run(["samtools", "faidx", fa_path], check=True)
    return fa_path


def check_dependencies(build_index, download_reference):
    """Verify required executables are on PATH."""
    required = ["generateMap.pl", "mapCounter", "bowtie"]
    if build_index:
        required.append("bowtie-build")
    if download_reference:
        required.append("samtools")
    missing = [cmd for cmd in required if shutil.which(cmd) is None]
    if missing:
        sys.exit(
            f"Error: the following commands were not found on PATH: {', '.join(missing)}\n"
            "Install hmmcopy_utils (https://github.com/shahcompbio/hmmcopy_utils) "
            "and bowtie (https://bowtie-bio.sourceforge.net)."
        )


def run_bowtie_build(reference, index_base, threads):
    """Build a bowtie index from the reference FASTA.

    Skips if the index files already exist.
    """
    if os.path.isfile(f"{index_base}.1.ebwt") or os.path.isfile(f"{index_base}.1.ebwtl"):
        print(f"[Step 0] [skip] Bowtie index already exists at {index_base}")
        return
    print(f"[Step 0] Building bowtie index at {index_base} ...")
    subprocess.run(
        ["bowtie-build", reference, index_base, "--threads", str(threads)],
        check=True,
    )
    print("  Done.")


def run_generate_map(reference, index_base, map_bw, read_length):
    """Run generateMap.pl to produce a per-base mappability BigWig.

    Skips if the output BigWig already exists.
    """
    if os.path.isfile(map_bw):
        print(f"[Step 1] [skip] Mappability BigWig already exists: {map_bw}")
        return
    print(f"[Step 1] Running generateMap.pl (read_length={read_length}) ...")
    result = subprocess.run(
        [
            "generateMap.pl",
            "-b", index_base,
            reference,
            "-l", str(read_length),
            "-o", map_bw,
        ],
        capture_output=True,
        text=True,
    )
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")
    if result.returncode != 0:
        sys.exit(f"generateMap.pl failed (exit code {result.returncode})")
    if not os.path.isfile(map_bw):
        sys.exit(
            f"generateMap.pl exited successfully but output file is missing: {map_bw}\n"
            "Check the stderr output above for bowtie errors."
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
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"mapCounter failed (exit code {result.returncode}):", file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr, end="")
        sys.exit(1)

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
            "Examples:\n"
            "  # With a local reference FASTA\n"
            "  python build_blacklist.hg38.py \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --out_file /path/to/blacklist.hg38.bed.gz \\\n"
            "    --build_index --threads 8\n\n"
            "  # Auto-download NCBI GRCh38 (downloaded to temp dir, removed on exit)\n"
            "  python build_blacklist.hg38.py \\\n"
            "    --download_reference \\\n"
            "    --out_file /path/to/blacklist.hg38.bed.gz \\\n"
            "    --threads 8"
        ),
    )
    ref_group = parser.add_mutually_exclusive_group(required=True)
    ref_group.add_argument(
        "--reference",
        help="Path to reference FASTA (must be indexed: .fai).",
    )
    ref_group.add_argument(
        "--download_reference", action="store_true", default=False,
        help="Download NCBI GRCh38 analysis set (no ALT, UCSC chr-names) to work dir. "
        "Removed on exit unless --work_dir is specified.",
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

    # Validate reference (if provided)
    if args.reference and not os.path.isfile(args.reference):
        sys.exit(f"Error: reference file not found: {args.reference}")

    # Check dependencies (download implies build_index)
    check_dependencies(args.build_index or args.download_reference, args.download_reference)

    # Set up work directory
    user_work_dir = args.work_dir is not None
    if user_work_dir:
        work_dir = os.path.abspath(args.work_dir)
        os.makedirs(work_dir, exist_ok=True)
    else:
        work_dir = tempfile.mkdtemp(prefix="blacklist_")

    # Resolve reference
    if args.download_reference:
        args.build_index = True  # downloaded FASTA has no pre-built index
        reference = download_reference(work_dir)
    else:
        reference = os.path.abspath(args.reference)

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
    index_base = os.path.join(work_dir, "ref_index")
    map_bw = os.path.join(work_dir, "mappability.bw")
    binned_tsv = os.path.join(work_dir, "binned_mappability.tsv")

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
    run_generate_map(reference, index_base, map_bw, args.read_length)

    # Step 2: mapCounter → binned DataFrame (cached to TSV)
    if os.path.isfile(binned_tsv):
        print(f"[Step 2] [skip] Loading cached binned mappability: {binned_tsv}")
        binned = pd.read_csv(binned_tsv, sep="\t")
    else:
        binned = run_map_counter(map_bw, args.window, STANDARD_CHROMS)
        binned.to_csv(binned_tsv, sep="\t", index=False)
        print(f"  Cached binned mappability to {binned_tsv}")

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

    # Step 5: merge and write
    combined = pd.concat(beds, ignore_index=True)
    merged = pr.PyRanges(combined).merge().sort()

    merged_df = merged.as_df()
    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    merged_df[["Chromosome", "Start", "End"]].to_csv(
        args.out_file, sep="\t", header=False, index=False,
        compression="gzip" if args.out_file.endswith(".gz") else None,
    )
    total_masked_bp = (merged_df["End"] - merged_df["Start"]).sum()
    print(
        f"[Step 5] Wrote {len(merged_df)} intervals to {args.out_file} "
        f"({total_masked_bp / 1e3:.1f} kbp masked)"
    )

    # Clean up auto-created temp dir only on success
    if not user_work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
