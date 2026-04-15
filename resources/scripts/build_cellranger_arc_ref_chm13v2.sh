#!/usr/bin/env bash
set -euo pipefail

# Build a Cell Ranger ARC reference for CHM13v2.0 (T2T) using UCSC hs1
# (chr*-style contig names).
#
# Based on: https://kb.10xgenomics.com/s/article/29207065679501-Building-a-Custom-T2T-reference-for-Cell-Ranger-ARC
#
# Requires: cellranger-arc (mkgtf, mkref), wget, gzip

FASTA_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
GTF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz"
MOTIFS_URL="https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <output_dir> [memgb]" >&2
  echo "  output_dir: directory to build the reference in" >&2
  echo "  memgb: memory in GB for mkref (default: 16)" >&2
  exit 1
fi

OUT=$1
MEMGB=${2:-16}
mkdir -p "${OUT}"
cd "${OUT}"

##################################################
# Step 1: Download FASTA, GTF, and motifs
##################################################
echo "=== Step 1: Download ==="

if [ ! -f "hs1.fa.gz" ]; then
  wget -c "${FASTA_URL}"
fi
if [ ! -f "hs1.ncbiRefSeq.gtf.gz" ]; then
  wget -c "${GTF_URL}"
fi
if [ ! -f "motifs.txt" ]; then
  wget -c -O motifs.txt "${MOTIFS_URL}"
fi

##################################################
# Step 2: Decompress
##################################################
echo "=== Step 2: Decompress ==="

if [ ! -f "hs1.fa" ]; then
  gunzip -k hs1.fa.gz
fi
if [ ! -f "hs1.ncbiRefSeq.gtf" ]; then
  gunzip -k hs1.ncbiRefSeq.gtf.gz
fi

##################################################
# Step 3: Filter annotations using mkgtf
##################################################
echo "=== Step 3: Filter GTF ==="

if [ ! -f "hs1.filtered.gtf" ]; then
  cellranger-arc mkgtf hs1.ncbiRefSeq.gtf hs1.filtered.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA
fi

##################################################
# Step 4: Generate reference using mkref
##################################################
echo "=== Step 4: Build reference ==="

cat > hs1.config <<'CONFIGEOF'
{
    organism: "human"
    genome: ["hs1"]
    input_fasta: ["hs1.fa"]
    input_gtf: ["hs1.filtered.gtf"]
    input_motifs: "motifs.txt"
}
CONFIGEOF

if [ ! -d "hs1" ]; then
  cellranger-arc mkref --config=hs1.config --memgb="${MEMGB}"
fi

echo "=== Done ==="
echo "Reference: ${OUT}/hs1"
echo "FASTA: ${OUT}/hs1/fasta/genome.fa"
echo "GTF: ${OUT}/hs1/genes/genes.gtf.gz"
echo "Use as gtf_file and reference in pipeline config."
