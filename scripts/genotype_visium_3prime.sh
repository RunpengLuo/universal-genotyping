#!/bin/bash
set -euo pipefail

# input
SAMPLE_FILE=$1
CONFIG_FILE=$2
OUTDIR=$3

source ${CONFIG_FILE}

########################################
TMPDIR=${OUTDIR}/tmp
LOGDIR=${OUTDIR}/log
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "$SCRIPT_DIR"

mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR}
mkdir -p ${LOGDIR}

mkdir -p ${LOGDIR}/genotype
mkdir -p ${LOGDIR}/allele_count
mkdir -p ${LOGDIR}/phase


SAMPLES=()
RANGER_DIRS=()
BAMS=()
while IFS=$'\t' read -r sample_name ranger_path; do
    echo "${sample_id}: ${ranger_path}"
    SAMPLES+=("$sample_id")
    RANGER_DIRS+=("$ranger_path")
done < ${SAMPLE_FILE}

# if [[ "$NO_NORMAL" == 0 ]]; then
#     [[ "${SAMPLES[0]}" == "normal" ]] || { echo "no normal sample found"; exit 1; }
#     NORMAL_BAM="${BAMS[0]}"
#     TUMOR_SAMPLES=("${SAMPLES[@]:1}")
#     TUMOR_BAMS=("${BAMS[@]:1}")
# fi
# NUM_TUMOR_SAMPLES=${#TUMOR_SAMPLES[@]}
NUM_SAMPLES=${#SAMPLES[@]}

# prepare bamfiles.txt and samples.txt
bamfile="${OUTDIR}/bams.txt"
samplefile="${OUTDIR}/samples.txt"
for ((i=0; i<${NUM_SAMPLES}; i++)); do
    echo "${SAMPLES[i]}" >> ${samplefile}
    echo "${RANGER_DIRS[i]}/possorted_genome_bam.bam" >> ${bamfile}
done

########################################
echo "genotyping on pseudobulk"
date
snp_dir="${OUTDIR}/genotype"
mkdir -p ${snp_dir}
if [[ ! -f "${snp_dir}/cellSNP.base.vcf.gz" ]]; then
    cellsnp-lite \
        -S ${bamfile} \
        -O ${snp_dir} \
        -R ${DB_SNP} \
        -p ${MAXJOBS} \
        --minMAF ${minMAF} \
        --minCOUNT ${minCOUNT} \
        --UMItag Auto \
        --cellTAG None \
        --gzip >> &>"${LOGDIR}/genotype/cellsnp_lite.log" &
else
    echo "skip"
fi

########################################
echo "label hom/het SNPs"
snp_vcf_file=${snp_dir}/snps.vcf.gz
if [[ ! -f ${snp_vcf_file} ]]; then
    python ${SCRIPT_DIR}/genotype_snps.py \
        -c ${snp_dir} \
        -o ${snp_vcf_file}
    tabix -f 
else

fi



# 1. panel SNPs VCF


        shapeit_concat_file="${TMPDIR}/shapeit.concat.lst"
        > ${shapeit_concat_file}
        while read -r CHROM START END; do
            region="${CHROM}:${START}-${END}"
            if [[ ! " ${CHROMS[*]} " =~ " ${CHROM#chr} " ]]; then
                echo "Skipping ${CHROM}:${START}-${END}"
                continue
            fi
            count=$(bcftools view -H -r "${region}" ${het_snp_file} | wc -l)
            if (( count < 100 )); then
                echo "#SNPs=${count} are too low at ${region}, skipping"
                continue
            fi

            echo "Phasing ${region}, #SNPs=${count} ..."
            phase_common_static \
                --input ${het_snp_file} \
                --map $(GET_PANEL_GMAP "${CHROM}") \
                --reference $(GET_PANEL_BCF "${CHROM}") \
                --region ${region} \
                --thread ${numThreads} \
                --log "${LOGDIR}/phase/shapeit.${region}.log" \
                --output "${TMPDIR}/phased.${region}.bcf" &>/dev/null
            
            if [[ -f "${TMPDIR}/phased.${region}.bcf" ]]; then
                bcftools view -Ov "${TMPDIR}/phased.${region}.bcf" | bgzip > "${TMPDIR}/phased.${region}.vcf.gz"
                bcftools index -f "${TMPDIR}/phased.${region}.vcf.gz"
                echo "${TMPDIR}/phased.${region}.vcf.gz" >> ${shapeit_concat_file}
            else
                echo "failed to phase ${region} likely due to no SNPs, check logs"
            fi
        done < ${REGION_BED}
        bcftools concat --file-list ${shapeit_concat_file} -Oz -o ${phase_file}
