# rule pileup_snps_cellsnp_lite_cell:
#     input:
#         barcode=lambda wc: get_data[(wc.data_type, wc.rep_id)][0],
#         bam=lambda wc: get_data[(wc.data_type, wc.rep_id)][1],
#         ranger=lambda wc: get_data[(wc.data_type, wc.rep_id)][2],
#         snp_file=lambda wc: branch(
#             run_genotype_snps,
#             then="phase/phased_snps.vcf.gz",
#             otherwise=config["ref_snp_file"],
#         ),
#     output:
#         cellsnp_file="pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz",
#         sample_file="pileup/{data_type}_{rep_id}/cellSNP.samples.tsv",
#         tot_mat="pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx",
#         ad_mat="pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx",
#     wildcard_constraints:
#         data_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
#     threads: config["threads"]["pileup"]
#     params:
#         cellsnp_lite=config["cellsnp_lite"],
#         UMItag=lambda wc: branch(
#             wc.data_type == "scATAC",
#             then="None",
#             otherwise=config["params_cellsnp_lite"]["UMItag"],
#         ),
#         cellTAG=config["params_cellsnp_lite"]["cellTAG"],
#         minMAF=config["params_cellsnp_lite"]["minMAF"],
#         minCOUNT=config["params_cellsnp_lite"]["minCOUNT_pileup"],
#         bcftools=config["bcftools"],
#     log:
#         "logs/pileup_snps.{data_type}_{rep_id}.log",
#     shell:
#         r"""
#         cellsnp_dir=$(dirname "{output.cellsnp_file}")
#         {params.cellsnp_lite} \
#             -b "{input.barcode}" \
#             -s "{input.bam}" \
#             -O ${{cellsnp_dir}} \
#             -R "{input.snp_file}" \
#             -p {threads} \
#             --minMAF {params.minMAF} \
#             --minCOUNT {params.minCOUNT} \
#             --UMItag {params.UMItag} \
#             --cellTAG {params.cellTAG} \
#             --gzip > {log} 2>&1
#         """


rule pileup_snps_cellsnp_lite_bulk:
    input:
        bam=lambda wc: get_data[(wc.data_type, wc.rep_id)][1],
        snp_file=lambda wc: branch(
            run_genotype_snps,
            then="phase/phased_snps.vcf.gz",
            otherwise=config["ref_snp_file"],
        ),
    output:
        cellsnp_file="pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz",
        sample_file="pileup/{data_type}_{rep_id}/cellSNP.samples.tsv",
        tot_mat="pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx",
        ad_mat="pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx",
    wildcard_constraints:
        data_type="bulkDNA",
    threads: config["threads"]["pileup"]
    params:
        cellsnp_lite=config["cellsnp_lite"],
        UMItag="None",
        cellTAG="None",
        minMAF=config["params_cellsnp_lite"]["minMAF"],
        minCOUNT=config["params_cellsnp_lite"]["minCOUNT_pileup"],
    log:
        "logs/pileup_snps.{data_type}_{rep_id}.log",
    shell:
        r"""
        cellsnp_dir=$(dirname "{output.cellsnp_file}")
        {params.cellsnp_lite} \
            -s "{input.bam}" \
            -O ${{cellsnp_dir}} \
            -R "{input.snp_file}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag {params.UMItag} \
            --cellTAG {params.cellTAG} \
            --gzip > {log} 2>&1
        """
