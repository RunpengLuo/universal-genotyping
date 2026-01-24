rule pileup_snps_cellsnp_lite_cell:
    input:
        barcode=lambda wc: get_data[(wc.data_type, wc.rep_id)][0],
        bam=lambda wc: get_data[(wc.data_type, wc.rep_id)][1],
        ranger=lambda wc: get_data[(wc.data_type, wc.rep_id)][2],
        snp_file=lambda wc: ("phase/phased_het_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
    output:
        cellsnp_file="pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz",
    wildcard_constraints:
        data_type="^(sc_GEX|sc_ATAC|VISIUM|VISIUM_3PRIME)$",
    threads: config["threads"]["pileup"]
    params:
        cellsnp_lite=config["cellsnp_lite"],
        UMItag=lambda wc: ("None" if wc.data_type == "sc_ATAC" else config["params_cellsnp_lite"]["UMItag"]),
        cellTAG=config["params_cellsnp_lite"]["cellTAG"],
        minMAF=0.0,
        minCOUNT=0,
        bcftools=config["bcftools"],
    log:
        "logs/pileup_snps.{data_type}_{rep_id}.log"
    shell:
        r"""
        set -euo pipefail
        cellsnp_dir=$(dirname "{output.cellsnp_file}")
        {params.cellsnp_lite} \
            -b "{input.barcode}" \
            -s "{input.bam}" \
            -O ${cellsnp_dir} \
            -R "{input.snp_file}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag {params.UMItag} \
            --cellTAG {params.cellTAG} \
            --gzip >> "{log}" 2>&1
        """

rule pileup_snps_cellsnp_lite_bulk:
    input:
        bam=lambda wc: get_data[(wc.data_type, wc.rep_id)][1],
        snp_file=lambda wc: ("phase/phased_het_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
    output:
        cellsnp_file="pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz",
    wildcard_constraints:
        data_type="^(bulk_DNA)$"
    threads: config["threads"]["pileup"]
    params:
        cellsnp_lite=config["cellsnp_lite"],
        UMItag="None",
        cellTAG="None",
        minMAF=0,
        minCOUNT=0,
    log:
        "logs/pileup_snps.{data_type}_{rep_id}.log"
    shell:
        r"""
        set -euo pipefail
        cellsnp_dir=$(dirname "{output.cellsnp_file}")
        {params.cellsnp_lite} \
            -s "{input.bam}" \
            -O ${cellsnp_dir} \
            -R "{input.snp_file}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag {params.UMItag} \
            --cellTAG {params.cellTAG} \
            --gzip >> "{log}" 2>&1
        """

# rule pileup_snps_bcftools_bulk:
#     input:
#         bam=lambda wc: get_data[(wc.data_type, wc.rep_id)][1],
#         snp_file=lambda wc: ("phase/phased_het_snps.vcf.gz" if run_genotype_snps else ref_snp_file),
#     output:
#         vcf_file="pileup/{data_type}_{rep_id}/allele_counts.vcf.gz",
#     shell:
#         r"""
#         echo "TODO"
#         """
