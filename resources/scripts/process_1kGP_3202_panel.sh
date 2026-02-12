#!/usr/bin/env bash
set -euo pipefail

BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

OUT=$1
mkdir -p "${OUT}"/{raw,snps,target_positions,phasing_panel}

: > "${OUT}/snp_vcfs.lst"
for chr in $(seq 1 22) X; do
  date
  echo "chr$chr"
  if [[ "${chr}" == "X" ]]; then
    vcfgz="1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
  else
    vcfgz="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  fi
  vcfgz_tbi=${vcfgz}.tbi
  if [ ! -f "${OUT}/raw/${vcfgz}" ]; then
    echo "download ${BASE}/${vcfgz}"
    curl -fL -C - -o "${OUT}/raw/${vcfgz}" "${BASE}/${vcfgz}"
    curl -fL -C - -o "${OUT}/raw/${vcfgz_tbi}" "${BASE}/${vcfgz_tbi}"
  fi

  if [ ! -f "${OUT}/snps/chr${chr}.vcf.gz" ]; then
    date
    echo "extract SNP sites"
    bcftools view -v snps -G -Oz -o "${OUT}/snps/chr${chr}.vcf.gz" "${OUT}/raw/${vcfgz}"
    tabix -p vcf "${OUT}/snps/chr${chr}.vcf.gz"
  fi
  echo "${OUT}/snps/chr${chr}.vcf.gz" >> "${OUT}/snp_vcfs.lst"

  if [ ! -f "${OUT}/target_positions/target.chr${chr}.pos.gz" ]; then
    date
    echo "extract SNP target positions"
    # add -m2 -M2 to filter multi-allelic positions.
    bcftools query -f '%CHROM\t%POS\n' \
      "${OUT}/snps/chr${chr}.vcf.gz" | bgzip -c > "${OUT}/target_positions/target.chr${chr}.pos.gz"
    tabix -s1 -b2 -e2 "${OUT}/target_positions/target.chr${chr}.pos.gz"
  fi
  
  if [ ! -f "${OUT}/phasing_panel/chr${chr}.genotypes.bcf" ]; then
    date
    echo "convert to phasing panel BCF file"
    bcftools view -Ob -o "${OUT}/phasing_panel/chr${chr}.genotypes.bcf" \
      "${OUT}/raw/${vcfgz}"
    bcftools index --csi "${OUT}/phasing_panel/chr${chr}.genotypes.bcf"
  fi
done

bcftools concat -f "${OUT}/snp_vcfs.lst" -Oz -o "${OUT}/snps.vcf.gz"
tabix -p vcf "${OUT}/snps.vcf.gz"

date
echo "=====summary====="
echo "raw files: ${OUT}/raw/*.vcf.gz"
echo "per-chromosome snp panel: ${OUT}/snps/chr<*>.vcf.gz"
echo "per-chromosome snp target positions: ${OUT}/target_positions/target.chr<*>.pos.gz"
echo "per-chromosome phasing panel: ${OUT}/phasing_panel/chr<*>.genotypes.bcf"
echo "merged snp panel: ${OUT}/snps.vcf.gz"
exit 0
