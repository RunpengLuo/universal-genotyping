#!/usr/bin/env python3
"""Build a per-window bias correction BED for WES data.

Reads a WES capture targets BED, merges overlapping intervals, adaptively
tiles each target into windows of at least ``--window`` bp, then computes
per-window GC content with optional mappability and replication timing columns.

Pipeline:
  1. Read & merge overlapping targets, filter to standard chroms
  2. Tile each merged target into windows via _tile_region()
  3. Assign gene name from parent target to each tiled window
  4. Filter blacklist + outside regions (subtract_blacklist)
  5. Assign region_id, compute GC/MAP/REPLI
  6. Sort and write output

Output is a gzipped TSV:
  #CHR  START  END  region_id  GC  [MAP]  [REPLI]  gene

Dependencies:
  - pybedtools, pandas, numpy
  - bigWigToBedGraph, liftOver (UCSC tools) -- only if --repliseq
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

from window_bed_utils import (
    REPLISEQ_REFVERS,
    _sort_windows,
    _tile_region,
    assign_region_id,
    compute_gc,
    compute_mappability,
    compute_repliseq,
    get_standard_chroms,
    print_summary,
    subtract_blacklist,
)


# ---------------------------------------------------------------------------
# Window generation (WES via adaptive tiling)
# ---------------------------------------------------------------------------


def generate_wes_windows(wes_targets_bed, window_size, standard_chroms, blacklist_bed):
    """Generate target windows for WES by merging and tiling capture targets.

    1. Read WES targets BED (BED4: chr/start/end/gene).
    2. Merge overlapping intervals, concatenating gene names.
    3. Filter to standard chromosomes.
    4. Tile each merged target into windows of at least ``window_size`` bp.
    5. Subtract blacklist regions.
    """
    # Read targets BED — support BED3 or BED4+
    targets = pd.read_csv(
        wes_targets_bed,
        sep="\t",
        header=None,
        comment="#",
    )
    if targets.shape[1] >= 4:
        targets = targets.iloc[:, :4]
        targets.columns = ["Chromosome", "Start", "End", "Name"]
    else:
        targets.columns = ["Chromosome", "Start", "End"]
        targets["Name"] = "."

    targets = targets[targets["Chromosome"].isin(standard_chroms)].reset_index(
        drop=True
    )
    n_raw = len(targets)
    raw_sizes = targets["End"] - targets["Start"]
    print(f"  {n_raw} raw target intervals on standard chromosomes")
    print(
        f"  raw target size: min={raw_sizes.min()}, mean={raw_sizes.mean():.0f}, "
        f"median={int(raw_sizes.median())}, max={raw_sizes.max()}"
    )

    # Merge overlapping intervals and collect gene names per merged interval
    targets = targets.sort_values(["Chromosome", "Start"]).reset_index(drop=True)
    merged_rows = []
    for chrom, grp in targets.groupby("Chromosome", sort=False):
        starts = grp["Start"].values
        ends = grp["End"].values
        names = grp["Name"].values
        cur_start, cur_end = starts[0], ends[0]
        cur_names = [names[0]]
        for i in range(1, len(starts)):
            if starts[i] <= cur_end:
                cur_end = max(cur_end, ends[i])
                cur_names.append(names[i])
            else:
                gene_str = (
                    ",".join(
                        sorted(set(n for n in cur_names if n != "." and pd.notna(n)))
                    )
                    or "."
                )
                merged_rows.append([chrom, cur_start, cur_end, gene_str])
                cur_start, cur_end = starts[i], ends[i]
                cur_names = [names[i]]
        gene_str = (
            ",".join(sorted(set(n for n in cur_names if n != "." and pd.notna(n))))
            or "."
        )
        merged_rows.append([chrom, cur_start, cur_end, gene_str])
    merged_df = pd.DataFrame(
        merged_rows, columns=["Chromosome", "Start", "End", "gene"]
    )

    merged_sizes = merged_df["End"] - merged_df["Start"]
    print(f"  {len(merged_df)} merged target intervals")
    print(
        f"  merged target size: min={merged_sizes.min()}, mean={merged_sizes.mean():.0f}, "
        f"median={int(merged_sizes.median())}, max={merged_sizes.max()}"
    )

    # Tile each merged target
    rows = []
    gene_labels = []
    chroms = merged_df["Chromosome"].values
    starts = merged_df["Start"].values.astype(int)
    ends = merged_df["End"].values.astype(int)
    genes = merged_df["gene"].values
    for i in range(len(merged_df)):
        tiled = _tile_region(chroms[i], starts[i], ends[i], window_size)
        rows.extend(tiled)
        gene_labels.extend([genes[i]] * len(tiled))

    windows = pd.DataFrame(rows, columns=["#CHR", "START", "END"])
    windows["gene"] = gene_labels

    # Subtract blacklist
    if blacklist_bed:
        n_before = len(windows)
        windows = subtract_blacklist(windows, blacklist_bed)
        n_removed = n_before - len(windows)
        print(
            f"  {len(windows)} windows after blacklist subtraction (removed {n_removed})"
        )

    sizes = windows["END"] - windows["START"]
    print(f"  {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes")
    print(
        f"  window size: min={sizes.min()}, mean={sizes.mean():.0f}, "
        f"median={int(sizes.median())}, max={sizes.max()}"
    )
    return windows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a per-window bias correction BED for WES. "
            "Merges and tiles WES capture targets into windows, then computes "
            "per-window GC content with optional mappability and replication timing."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python build_wes_window_bed.py \\\n"
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
        help="Output gzipped TSV path (e.g. gc_map.wes.bed.gz).",
    )
    # Optional filtering
    parser.add_argument(
        "--blacklist_bed",
        default=None,
        help="Optional blacklist BED -- subtract these regions.",
    )
    # WES options
    parser.add_argument(
        "--wes_targets_bed",
        required=True,
        help="Vendor exon capture targets BED (BED3 or BED4 with gene name).",
    )
    parser.add_argument(
        "--window", type=int, default=267, help="Window size in bp (default: 267)."
    )
    # Covariate options
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
        help="hg19-to-hg38 liftOver chain file (downloaded if omitted with --repliseq).",
    )
    parser.add_argument(
        "--bigwig_dir",
        default=None,
        help="Directory for pre-downloaded WaveSignal bigWig files.",
    )
    parser.add_argument(
        "--work_dir",
        default=None,
        help="Working directory for intermediate files (default: temp dir).",
    )
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

    t0 = time.time()

    print("=== build_wes_window_bed.py ===")
    print(f"  Ref version: {args.reference_version}")
    print(f"  Reference:   {args.reference}")
    print(f"  Genome size: {args.genome_size}")
    print(f"  Region BED:  {args.region_bed}")
    print(f"  Blacklist:   {args.blacklist_bed or 'none'}")
    print(f"  WES targets: {args.wes_targets_bed}")
    print(f"  Window:      {args.window} bp")
    print(f"  Mappability: {args.mappability_bed or 'none'}")
    print(f"  Repli-seq:   {'yes' if args.repliseq else 'no'}")
    print(f"  Output:      {args.out_file}")
    print()

    standard_chroms = get_standard_chroms(args.reference_version, args.genome_size)
    print(f"  {len(standard_chroms)} standard chromosomes")

    # --- Step 1: Generate windows via adaptive tiling ---
    print("[1/6] WES: merging targets and tiling windows ...")
    windows = generate_wes_windows(
        args.wes_targets_bed,
        args.window,
        standard_chroms,
        args.blacklist_bed,
    )

    # --- Step 2: Assign region_id, drop windows outside regions ---
    print("[2/6] Assigning region_id ...")
    windows["region_id"] = assign_region_id(windows, args.region_bed)
    n_before = len(windows)
    unassigned_mask = windows["region_id"].isna()
    n_dropped = int(unassigned_mask.sum())
    if n_dropped > 0:
        print(
            f"  WARNING: {n_dropped}/{n_before} windows could not be assigned a region_id"
        )
        print(windows.loc[unassigned_mask].head(30).to_string(index=False))
        windows = windows[~unassigned_mask].reset_index(drop=True)
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
            windows,
            args.mappability_bed,
            args.genome_size,
        )
        out_cols.append("MAP")
    else:
        print("[4/6] Skipping mappability (no --mappability_bed)")

    # --- Step 5: Compute replication timing (optional) ---
    if args.repliseq:
        if args.reference_version not in REPLISEQ_REFVERS:
            print(
                f"[5/6] Skipping repli-seq (unsupported for {args.reference_version})"
            )
        else:
            print("[5/6] Computing replication timing ...")
            windows["REPLI"] = compute_repliseq(windows, standard_chroms, args)
            out_cols.append("REPLI")
    else:
        print("[5/6] Skipping repli-seq (--repliseq not set)")

    out_cols.append("gene")

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
