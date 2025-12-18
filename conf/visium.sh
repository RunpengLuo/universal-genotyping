# general parameters
READ_TYPE=VISIUM # TGS,NGS,WES
NO_NORMAL=0 # TODO relax later
MAXJOBS=4

# cellsnp-lite parameters
minMAF=0
minCOUNT=2



# # bcftools pileup parameters
# q=0
# Q=11
# d=300
# mincov=5
# CHROMS=($(seq 1 22))

# # # SNP filtering parameters
# # gamma=0.05
# # min_ad=1

# # # mosdpeth parameters
# # readquality=11

# # phasing parameters
# minMAPQ=20
# numThreads=8

# chm13v2 hg38
reference_version=chm13v2
if [ "${reference_version}" == "chm13v2" ]; then
    # T2T-CHM13v2.0
    REGION_BED=/diskmnt/Users2/runpengl/data/chm13v2.0_region.bed
    DB_SNP=/diskmnt/Projects/ccRCC_longread/runpengl/vcf/chm13v2.0_dbSNPv155.vcf.gz
    REFERENCE=/diskmnt/Projects/ccRCC_longread/runpengl/reference/T2T-CHM13v2.0.fasta
    GENOME_SIZE=/diskmnt/Projects/ccRCC_longread/runpengl/reference/T2T-CHM13v2.0.sizes
    REF_PANEL=/diskmnt/Projects/ccRCC_longread/runpengl/data/1000G_chm13v2/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf.gz
    GEN_MAP_DIR=/diskmnt/Projects/ccRCC_longread/runpengl/data/1000G_chm13v2/phasing_T2T/resources/recombination_maps/t2t_native_scaled_maps
elif [ "${reference_version}" == "hg38" ]; then
    # GRCh38
    REGION_BED=/diskmnt/Projects/ccRCC_longread/runpengl/reference/hg38_region.bed
    DB_SNP=/diskmnt/Projects/ccRCC_longread/runpengl/data/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.chr_sorted.vcf.gz
    REFERENCE=/diskmnt/Projects/ccRCC_longread/runpengl/reference/GRCh38.fasta
    GENOME_SIZE=/diskmnt/Projects/ccRCC_longread/runpengl/reference/GRCh38.sizes
    REF_PANEL=/diskmnt/Projects/ccRCC_longread/runpengl/data/1000G_hg38
    GEN_MAP_DIR=/diskmnt/Projects/ccRCC_longread/runpengl/data/shapeit5_resources/maps/b38
fi

GET_PANEL_BCF() {
	local chrom="$1"
    if [ "${reference_version}" == "chm13v2" ]; then
        panel_bcf="${REF_PANEL}"
    elif [ "${reference_version}" == "hg38" ]; then
        panel_bcf="${REF_PANEL}/${CHROM}.genotypes.bcf"
    fi
    echo ${panel_bcf}
}

GET_PANEL_GMAP() {
	local chrom="$1"
    if [ "${reference_version}" == "chm13v2" ]; then
        panel_gmap="${GEN_MAP_DIR}/${CHROM}.t2t.scaled.gmap.gz"
    elif [ "${reference_version}" == "hg38" ]; then
        panel_gmap="${GEN_MAP_DIR}/${CHROM}.b38.gmap.gz"
    fi
    echo ${panel_gmap}
}
