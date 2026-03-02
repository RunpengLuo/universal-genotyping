#!/usr/bin/env python3
"""Download UCSC hg38 repeat annotations and low-mappability regions, merge into a single blacklist BED."""

import argparse
import os
import sys
import tempfile
import urllib.request

import numpy as np
import pandas as pd
import pyBigWig
import pyranges as pr

UCSC_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"
TABLES = ["genomicSuperDups.txt.gz", "simpleRepeat.txt.gz"]
MAPPABILITY_URL = (
    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/"
    "k100.Umap.MultiTrackMappability.bw"
)
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


def get_low_mappability_regions(bw_path, min_map_score, chroms):
    """Extract regions with mappability score < min_map_score from a bigWig file.

    Regions not covered by the bigWig (score=0) are also included.
    """
    bw = pyBigWig.open(bw_path)
    rows = []
    sorted_chroms = sorted(chroms)
    for i, chrom in enumerate(sorted_chroms):
        print(f"  Processing {chrom} ({i+1}/{len(sorted_chroms)})");
        chrom_len = bw.chroms().get(chrom)
        if chrom_len is None:
            continue
        intervals = bw.intervals(chrom)
        if intervals is None:
            rows.append((chrom, 0, chrom_len))
            continue
        arr = np.array(intervals)
        starts = arr[:, 0].astype(np.int64)
        ends = arr[:, 1].astype(np.int64)
        scores = arr[:, 2]

        # low-score intervals
        low = scores < min_map_score
        if low.any():
            chrom_col = np.full(low.sum(), chrom)
            rows.append(np.column_stack([chrom_col, starts[low], ends[low]]))

        # gaps between consecutive intervals (unmappable, score=0)
        gap_starts = ends[:-1]
        gap_ends = starts[1:]
        gap_mask = gap_ends > gap_starts
        if gap_mask.any():
            chrom_col = np.full(gap_mask.sum(), chrom)
            rows.append(np.column_stack([chrom_col, gap_starts[gap_mask], gap_ends[gap_mask]]))

        # leading/trailing gaps
        if starts[0] > 0:
            rows.append(np.array([[chrom, 0, starts[0]]]))
        if ends[-1] < chrom_len:
            rows.append(np.array([[chrom, ends[-1], chrom_len]]))
    bw.close()
    if not rows:
        return pd.DataFrame(columns=["Chromosome", "Start", "End"])
    combined = np.vstack(rows)
    return pd.DataFrame(combined, columns=["Chromosome", "Start", "End"]).astype(
        {"Start": np.int64, "End": np.int64}
    )


def main():
    parser = argparse.ArgumentParser(
        description="Build a merged repeat/segdup/low-mappability blacklist BED for hg38."
    )
    parser.add_argument(
        "--out_file", required=True, help="Output BED file path."
    )
    parser.add_argument(
        "--min_map_score",
        type=float,
        default=0.9,
        help="Exclude regions with mappability score below this threshold (default: 0.9).",
    )
    args = parser.parse_args()

    beds = []
    tmpdir = tempfile.mkdtemp()
    tmp_files = []
    try:
        # UCSC repeat/segdup tables
        for table in TABLES:
            local = download_file(f"{UCSC_BASE}/{table}", tmpdir, table)
            tmp_files.append(local)
            bed = read_bed_columns(local)
            print(f"  {table}: {len(bed)} intervals")
            beds.append(bed)

        # UCSC mappability bigWig
        bw_path = download_file(MAPPABILITY_URL, tmpdir)
        tmp_files.append(bw_path)
        low_map = get_low_mappability_regions(bw_path, args.min_map_score, STANDARD_CHROMS)
        print(f"Low-mappability regions (score < {args.min_map_score}): {len(low_map)} intervals")
        beds.append(low_map)
    finally:
        for p in tmp_files:
            if os.path.exists(p):
                os.remove(p)
        os.rmdir(tmpdir)

    combined = pd.concat(beds, ignore_index=True)
    merged = pr.PyRanges(combined).merge().sort()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    merged.as_df()[["Chromosome", "Start", "End"]].to_csv(
        args.out_file, sep="\t", header=False, index=False
    )
    print(f"Wrote {len(merged)} intervals to {args.out_file}")


if __name__ == "__main__":
    main()
