#!/usr/bin/env bash
set -euo pipefail

# Prepare 1kGP n=3,202 phased panel for CHM13v2.0 (T2T) reference.
# Splits the whole-genome biallelic BCF into per-chromosome SNP panel VCFs,
# target position files, and phasing panel BCFs.
#
# If --bcf is provided, uses that file directly (no download).
# Otherwise downloads from HPRC S3.
#
# Requires: bcftools, bgzip, tabix

S3_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1"
BCF_NAME="1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf.gz"

usage() {
  echo "Usage: $0 <output_dir> [--bcf /path/to/existing.bcf.gz]" >&2
  exit 1
}

[[ $# -lt 1 ]] && usage

OUT=$1; shift
BCF_INPUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bcf) BCF_INPUT="$2"; shift 2 ;;
    *) usage ;;
  esac
done

mkdir -p "${OUT}"/{snps,target_positions,phasing_panel}

if [[ -n "${BCF_INPUT}" ]]; then
  if [[ ! -f "${BCF_INPUT}" ]]; then
    echo "ERROR: BCF file not found: ${BCF_INPUT}" >&2
    exit 1
  fi
  echo "using existing BCF: ${BCF_INPUT}"
  INPUT_BCF="${BCF_INPUT}"
else
  mkdir -p "${OUT}/raw"
  INPUT_BCF="${OUT}/raw/${BCF_NAME}"
  if [ ! -f "${INPUT_BCF}" ]; then
    date
    echo "downloading ${BCF_NAME}"
    wget -c -P "${OUT}/raw" "${S3_BASE}/${BCF_NAME}"
    wget -c -P "${OUT}/raw" "${S3_BASE}/${BCF_NAME}.csi"
  fi
fi

: > "${OUT}/snp_vcfs.lst"
for chr in $(seq 1 22) X; do
  date
  echo "chr${chr}"

  if [ ! -f "${OUT}/snps/chr${chr}.vcf.gz" ]; then
    echo "extract SNP sites for chr${chr}"
    bcftools view -v snps -G -r "chr${chr}" -Oz \
      -o "${OUT}/snps/chr${chr}.vcf.gz" "${INPUT_BCF}"
    tabix -p vcf "${OUT}/snps/chr${chr}.vcf.gz"
  fi
  echo "${OUT}/snps/chr${chr}.vcf.gz" >> "${OUT}/snp_vcfs.lst"

  if [ ! -f "${OUT}/target_positions/target.chr${chr}.pos.gz" ]; then
    echo "extract target positions for chr${chr}"
    bcftools query -f '%CHROM\t%POS\n' \
      "${OUT}/snps/chr${chr}.vcf.gz" | bgzip -c > "${OUT}/target_positions/target.chr${chr}.pos.gz"
    tabix -s1 -b2 -e2 "${OUT}/target_positions/target.chr${chr}.pos.gz"
  fi

  if [ ! -f "${OUT}/phasing_panel/chr${chr}.genotypes.bcf" ]; then
    echo "extract phasing panel BCF for chr${chr}"
    bcftools view -r "chr${chr}" -Ob \
      -o "${OUT}/phasing_panel/chr${chr}.genotypes.bcf" "${INPUT_BCF}"
    bcftools index --csi "${OUT}/phasing_panel/chr${chr}.genotypes.bcf"
  fi
done

bcftools concat -f "${OUT}/snp_vcfs.lst" -Oz -o "${OUT}/snps.vcf.gz"
tabix -p vcf "${OUT}/snps.vcf.gz"

date
echo "=====summary====="
echo "input BCF: ${INPUT_BCF}"
echo "per-chromosome snp panel: ${OUT}/snps/chr<*>.vcf.gz"
echo "per-chromosome snp target positions: ${OUT}/target_positions/target.chr<*>.pos.gz"
echo "per-chromosome phasing panel: ${OUT}/phasing_panel/chr<*>.genotypes.bcf"
echo "merged snp panel: ${OUT}/snps.vcf.gz"
exit 0
