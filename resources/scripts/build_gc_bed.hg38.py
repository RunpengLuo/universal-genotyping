#!/usr/bin/env python3
"""Build a fixed-window GC-content (+ optional mappability) BED for hg38.

Generates non-overlapping windows across a reference genome and computes
per-window GC fraction using pybedtools ``nucleotide_content()``.  If a
mappability BigWig/BED is provided, per-window mean mappability is also
computed via ``BedTool.map()``.

Output is a gzipped TSV with columns: #CHR, START, END, GC, MAP.

Dependencies:
  - pybedtools
  - pandas, numpy
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from pybedtools import BedTool


STANDARD_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


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
    map_vals = pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(1.0).clip(0.0, 1.0)
    return map_vals


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a fixed-window GC-content BED for a reference genome. "
            "Optionally includes per-window mean mappability from a BED track."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # GC-only (500 bp windows)\n"
            "  python build_gc_bed.hg38.py \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --out_file gc_500bp.hg38.bed.gz\n\n"
            "  # GC + mappability\n"
            "  python build_gc_bed.hg38.py \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --mappability_bed /path/to/mappability.bed \\\n"
            "    --out_file gc_map_500bp.hg38.bed.gz"
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
        help="Output gzipped TSV path (e.g. gc_500bp.hg38.bed.gz).",
    )
    parser.add_argument(
        "--window", type=int, default=500,
        help="Window size in bp (default: 500).",
    )
    parser.add_argument(
        "--mappability_bed", default=None,
        help="Optional BED-format mappability track (4th column = score).",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.reference):
        sys.exit(f"Error: reference not found: {args.reference}")
    if not os.path.isfile(args.genome_size):
        sys.exit(f"Error: genome_size not found: {args.genome_size}")
    if args.mappability_bed and not os.path.isfile(args.mappability_bed):
        sys.exit(f"Error: mappability_bed not found: {args.mappability_bed}")

    print("=== build_gc_bed.hg38.py ===")
    print(f"  Reference:   {args.reference}")
    print(f"  Genome size: {args.genome_size}")
    print(f"  Window:      {args.window} bp")
    print(f"  Mappability: {args.mappability_bed or 'none (MAP=1.0)'}")
    print(f"  Output:      {args.out_file}")
    print()

    # Step 1: generate windows
    print("[Step 1] Generating fixed windows ...")
    windows = generate_windows(args.genome_size, args.window)
    print(f"  {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes")

    # Step 2: compute GC
    print("[Step 2] Computing GC content ...")
    gc_df = compute_gc(windows, args.reference)
    print(f"  GC range: [{gc_df['GC'].min():.4f}, {gc_df['GC'].max():.4f}]")

    # Step 3: compute mappability (or fill 1.0)
    if args.mappability_bed:
        print("[Step 3] Computing mappability ...")
        gc_df["MAP"] = compute_mappability(gc_df, args.mappability_bed, args.genome_size)
    else:
        print("[Step 3] No mappability track; setting MAP=1.0")
        gc_df["MAP"] = 1.0

    # Step 4: write output
    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    gc_df[["#CHR", "START", "END", "GC", "MAP"]].to_csv(
        args.out_file, sep="\t", header=True, index=False,
        compression="gzip" if args.out_file.endswith(".gz") else None,
    )
    print(f"[Step 4] Wrote {len(gc_df)} rows to {args.out_file}")


if __name__ == "__main__":
    main()
