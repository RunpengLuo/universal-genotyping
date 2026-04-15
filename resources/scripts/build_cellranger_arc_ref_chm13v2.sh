#!/usr/bin/env bash
set -euo pipefail

# Build a Cell Ranger ARC reference for CHM13v2.0 (T2T).
# Following: https://kb.10xgenomics.com/s/article/29207065679501-Building-a-Custom-T2T-reference-for-Cell-Ranger-ARC
#
# Requires: cellranger-arc (mkgtf, mkref), wget, awk, gzip

NCBI_FTP="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0"
FASTA_GZ="GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"
GTF_GZ="GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz"
MOTIFS_URL="https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <output_dir> [memgb]" >&2
  echo "  output_dir: directory to build the reference in" >&2
  echo "  memgb: memory in GB for mkref (default: 96)" >&2
  exit 1
fi

OUT=$1
MEMGB=${2:-96}
mkdir -p "${OUT}"
cd "${OUT}"

##################################################
# Step 1: Download FASTA, GTF, and motifs
##################################################
echo "=== Step 1: Download ==="

if [ ! -f "${FASTA_GZ}" ]; then
  wget -c "${NCBI_FTP}/${FASTA_GZ}"
fi
if [ ! -f "${GTF_GZ}" ]; then
  wget -c "${NCBI_FTP}/${GTF_GZ}"
fi
if [ ! -f "motifs.txt" ]; then
  wget -c -O motifs.txt "${MOTIFS_URL}"
fi

# Decompress
if [ ! -f "t2t_genome.fa" ]; then
  echo "decompressing FASTA"
  gunzip -c "${FASTA_GZ}" > t2t_genome_raw.fa
fi
if [ ! -f "t2t_raw.gtf" ]; then
  echo "decompressing GTF"
  gunzip -c "${GTF_GZ}" > t2t_raw.gtf
fi

##################################################
# Step 2: Format FASTA
# - Retain sequence name only (remove descriptions)
# - Convert bases to upper case
##################################################
echo "=== Step 2: Format FASTA ==="

if [ ! -f "t2t_genome.fa" ]; then
  awk '{if ($1 ~ />/) {print $1} else {print $_}}' t2t_genome_raw.fa \
    | tr '[:lower:]' '[:upper:]' > t2t_genome.fa
  rm -f t2t_genome_raw.fa
fi

##################################################
# Step 3: Filter annotations using mkgtf
##################################################
echo "=== Step 3: Filter GTF ==="

if [ ! -f "t2t.gtf" ]; then
  cellranger-arc mkgtf t2t_raw.gtf t2t.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA
fi

##################################################
# Step 4: Generate reference using mkref
##################################################
echo "=== Step 4: Build reference ==="

cat > t2t.config <<'CONFIGEOF'
{
    organism: "human"
    genome: ["t2t"]
    input_fasta: ["t2t_genome.fa"]
    input_gtf: ["t2t.gtf"]
    input_motifs: "motifs.txt"
}
CONFIGEOF

if [ ! -d "t2t" ]; then
  cellranger-arc mkref --config=t2t.config --memgb="${MEMGB}"
fi

echo "=== Done ==="
echo "Reference: ${OUT}/t2t"
echo "FASTA: ${OUT}/t2t/fasta/genome.fa"
echo "GTF: ${OUT}/t2t/genes/genes.gtf.gz"
echo "Use as gtf_file in pipeline config."
