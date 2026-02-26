#!/usr/bin/env python3
"""Download UCSC hg38 repeat annotations and merge into a single blacklist BED."""

import argparse
import gzip
import os
import tempfile
import urllib.request

import pandas as pd
import pyranges as pr

UCSC_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database"
TABLES = ["genomicSuperDups.txt.gz", "simpleRepeat.txt.gz"]
STANDARD_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def download_table(name, dest_dir):
    """Download a UCSC table file to *dest_dir* and return the local path."""
    url = f"{UCSC_BASE}/{name}"
    local_path = os.path.join(dest_dir, name)
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, local_path)
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


def main():
    parser = argparse.ArgumentParser(
        description="Build a merged repeat/segdup blacklist BED for hg38."
    )
    parser.add_argument(
        "--out_file", required=True, help="Output BED file path."
    )
    args = parser.parse_args()

    beds = []
    tmpdir = tempfile.mkdtemp()
    try:
        for table in TABLES:
            local = download_table(table, tmpdir)
            beds.append(read_bed_columns(local))
    finally:
        for table in TABLES:
            p = os.path.join(tmpdir, table)
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
