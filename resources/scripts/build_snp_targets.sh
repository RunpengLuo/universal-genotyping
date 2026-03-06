#!/usr/bin/env bash
set -euo pipefail

# Build per-chromosome SNP target position files from any SNP panel VCF.
# Output: target.chr{1..22,X}.pos.gz + .tbi in the specified output directory.
#
# Usage: bash build_snp_targets.sh <input.vcf.gz> <output_dir>
#
# Requires: bcftools, bgzip, tabix

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.vcf.gz> <output_dir>" >&2
  exit 1
fi

INPUT_VCF=$1
OUT=$2

if [[ ! -f "${INPUT_VCF}" ]]; then
  echo "Error: input VCF not found: ${INPUT_VCF}" >&2
  exit 1
fi

mkdir -p "${OUT}"

for chr in $(seq 1 22) X; do
  outfile="${OUT}/target.chr${chr}.pos.gz"
  if [[ -f "${outfile}" && -f "${outfile}.tbi" ]]; then
    echo "skip chr${chr} (already exists)"
    continue
  fi
  echo "chr${chr}"
  bcftools query -r "chr${chr}" -f '%CHROM\t%POS\n' "${INPUT_VCF}" \
    | bgzip -c > "${outfile}"
  tabix -s1 -b2 -e2 "${outfile}"
done

echo "Done. Output: ${OUT}/target.chr{1..22,X}.pos.gz"
