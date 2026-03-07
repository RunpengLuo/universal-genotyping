#!/usr/bin/env python3
"""Build a per-window bias correction BED for WES data.

Wraps ``cnvkit.py target`` and ``cnvkit.py antitarget`` to generate
target + antitarget windows, then computes per-window GC content with
optional mappability and replication timing columns.

Pipeline:
  1. cnvkit.py target --split  (subdivide exon targets)
  2. cnvkit.py antitarget      (off-target bins in accessible regions)
  3. Combine, filter standard chroms, subtract blacklist
  4. Assign region_id, compute GC/MAP/REPLI
  5. Sort and write output

Output is a gzipped TSV:
  #CHR  START  END  region_id  GC  [MAP]  [REPLI]  is_target

Dependencies:
  - cnvkit.py (must be on PATH)
  - pybedtools, pyranges, pandas, numpy
  - bigWigToBedGraph, liftOver (UCSC tools) -- only if --repliseq
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

import pandas as pd

from window_bed_utils import (
    _sort_windows,
    assign_region_id,
    compute_gc,
    compute_mappability,
    compute_repliseq,
    get_standard_chroms,
    print_summary,
    subtract_blacklist,
)


# ---------------------------------------------------------------------------
# Window generation (WES via CNVkit)
# ---------------------------------------------------------------------------


def generate_wes_windows(wes_targets_bed, region_bed, blacklist_bed, standard_chroms):
    """Generate target + antitarget windows for WES using cnvkit.py."""
    if shutil.which("cnvkit.py") is None:
        sys.exit("Error: cnvkit.py not found on PATH")

    with tempfile.TemporaryDirectory() as tmpdir:
        # cnvkit.py target (--split, default avg_size=267)
        tgt_out = os.path.join(tmpdir, "targets.bed")
        subprocess.check_call([
            "cnvkit.py", "target", wes_targets_bed,
            "--split", "-o", tgt_out,
        ])

        # cnvkit.py antitarget (default avg=150kb, min=auto ~9374bp)
        anti_out = os.path.join(tmpdir, "antitargets.bed")
        cmd = ["cnvkit.py", "antitarget", tgt_out, "-o", anti_out]
        if region_bed:
            cmd.extend(["-g", region_bed])
        subprocess.check_call(cmd)

        # Read BED4 outputs, combine with is_target
        targets = pd.read_csv(tgt_out, sep="\t", header=None,
                              names=["#CHR", "START", "END", "gene"])
        antitargets = pd.read_csv(anti_out, sep="\t", header=None,
                                  names=["#CHR", "START", "END", "gene"])

    targets["is_target"] = 1
    antitargets["is_target"] = 0
    windows = pd.concat([targets, antitargets], ignore_index=True)
    windows = windows.drop(columns=["gene"])
    windows = windows[windows["#CHR"].isin(standard_chroms)].reset_index(drop=True)

    n_tgt = int(windows["is_target"].sum())
    n_anti = len(windows) - n_tgt
    print(f"  {len(windows)} windows ({n_tgt} target, {n_anti} antitarget)")
    return windows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a per-window bias correction BED for WES. "
            "Uses cnvkit.py target + antitarget to generate windows, then "
            "computes per-window GC content with optional mappability and "
            "replication timing."
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
    parser.add_argument("--reference_version", required=True,
                        help="Reference genome version (hg19, hg38, chm13v2).")
    parser.add_argument("--reference", required=True,
                        help="Path to reference FASTA (must have .fai index).")
    parser.add_argument("--genome_size", required=True,
                        help="Path to chromosome sizes file (tab-separated: chrom, size).")
    parser.add_argument("--region_bed", required=True,
                        help="Accessible regions BED -- keep only windows overlapping these regions.")
    parser.add_argument("--out_file", required=True,
                        help="Output gzipped TSV path (e.g. gc_map.wes.bed.gz).")
    # Optional filtering
    parser.add_argument("--blacklist_bed", default=None,
                        help="Optional blacklist BED -- subtract these regions.")
    # WES options
    parser.add_argument("--wes_targets_bed", required=True,
                        help="Vendor exon capture targets BED (required).")
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

    t0 = time.time()

    print("=== build_wes_window_bed.py ===")
    print(f"  Ref version: {args.reference_version}")
    print(f"  Reference:   {args.reference}")
    print(f"  Genome size: {args.genome_size}")
    print(f"  Region BED:  {args.region_bed}")
    print(f"  Blacklist:   {args.blacklist_bed or 'none'}")
    print(f"  WES targets: {args.wes_targets_bed}")
    print(f"  Mappability: {args.mappability_bed or 'none'}")
    print(f"  Repli-seq:   {'yes' if args.repliseq else 'no'}")
    print(f"  Output:      {args.out_file}")
    print()

    standard_chroms = get_standard_chroms(args.reference_version, args.genome_size)
    print(f"  {len(standard_chroms)} standard chromosomes")

    # --- Step 1: Generate windows via cnvkit.py ---
    print("[1/6] WES: generating target + antitarget windows via cnvkit.py ...")
    windows = generate_wes_windows(
        args.wes_targets_bed, args.region_bed, args.blacklist_bed, standard_chroms,
    )

    # Subtract blacklist from combined target+antitarget windows
    if args.blacklist_bed:
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
