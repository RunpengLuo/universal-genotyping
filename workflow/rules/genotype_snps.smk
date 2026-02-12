"""
Inputs
1. BAM files
2. SNP panels
3. reference genome

Outputs: bi-allelic hom-alt and het ref/alt SNPs, per chromosome.
snps/<chrom>.vcf.gz
"""

if workflow_mode == "bulk_genotyping":

    # TODO genotyping tumor samples
    # handle mulitple BAM files
    rule genotype_snps_bulk:
        input:
            bams=lambda wc: branch(
                len(normal_bams) > 0, then=normal_bams[0], otherwise=tumor_bams[0]
            ),
            target_pos=lambda wc: config["snp_targets"] + "/target.chr{chrname}.pos.gz",
            reference=config["reference"],
        output:
            snp_vcf=config["snp_dir"] + "/chr{chrname}.vcf.gz",
        log:
            config["log_dir"] + "/genotype_snps_bulk/chr{chrname}.log",
        threads: config["threads"]["genotype"]
        params:
            chrom="chr{chrname}",
            bcftools=config["bcftools"],
            min_mapq=config["params_bcftools"]["min_mapq"],
            min_baseq=config["params_bcftools"]["min_baseq"],
            min_dp=config["params_bcftools"]["min_dp"],
            max_depth=config["params_bcftools"]["max_depth"],
        shell:
            r"""
            {params.bcftools} mpileup {input.bams} \
                -f "{input.reference}" \
                -Ou \
                -a INFO/AD,AD,DP \
                --skip-indels \
                -q {params.min_mapq} \
                -Q {params.min_baseq} \
                -d {params.max_depth} \
                -T {input.target_pos} \
            | {params.bcftools} call -m -Ou \
            | {params.bcftools} view -v snps -m2 -M2 \
                -i 'GT="alt" && FMT/DP>={params.min_dp}' \
                -Oz -o {output.snp_vcf} > {log} 2>&1

            tabix -p vcf {output.snp_vcf}
            """


if workflow_mode == "single_cell_genotyping":

    rule genotype_snps_pseudobulk_mode1b:
        input:
            bams=lambda wc: modality2bams[wc.modality],
            snp_panel=config["snp_panel"],
        output:
            out_dir=directory(config["snp_dir"] + "/pseudobulk_{modality}"),
            out_vcf=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.base.vcf.gz",
            out_tsv=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.samples.tsv",
            out_dp=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.tag.DP.mtx",
            out_ad=config["snp_dir"] + "/pseudobulk_{modality}/cellSNP.tag.AD.mtx",
            bam_lst=temp("tmp/bams.{modality}.lst"),
        log:
            config["log_dir"] + "/genotype_snps_pseudobulk.{modality}.log",
        threads: config["threads"]["genotype"]
        params:
            bcftools=config["bcftools"],
            cellsnp_lite=config["cellsnp_lite"],
            UMItag=lambda wc: branch(
                wc.modality == "RNA",
                then=config["params_cellsnp_lite"]["UMItag"],
                otherwise="None",
            ),
            minMAF=config["params_cellsnp_lite"]["minMAF_genotype"],
            minCOUNT=config["params_cellsnp_lite"]["minCOUNT_genotype"],
        shell:
            r"""
            printf "%s\n" {input.bams} > "{output.bam_lst}"
            {params.cellsnp_lite} \
                -S "{output.bam_lst}" \
                -R "{input.snp_panel}" \
                -O "{output.out_dir}" \
                -p {threads} \
                --minMAF {params.minMAF} \
                --minCOUNT {params.minCOUNT} \
                --UMItag {params.UMItag} \
                --cellTAG None \
                --gzip > {log} 2>&1
            """

    rule annotate_snps_pseudobulk:
        input:
            raw_snp_vcfs=[
                config["snp_dir"] + f"/pseudobulk_{modality}/cellSNP.base.vcf.gz"
                for modality in modalities
            ],
            genome_size=config["genome_size"],
        output:
            snp_vcfs=expand(
                config["snp_dir"] + "/chr{chrname}.vcf.gz",
                chrname=config["chromosomes"],
            ),
            snp_vcfs_tbi=expand(
                config["snp_dir"] + "/chr{chrname}.vcf.gz.tbi",
                chrname=config["chromosomes"],
            ),
        params:
            modalities=modalities,
            min_het_reads=config["params_annotate_snps"]["min_het_reads"],
            min_hom_dp=config["params_annotate_snps"]["min_hom_dp"],
            min_vaf_thres=config["params_annotate_snps"]["min_vaf_thres"],
            filter_nz_OTH=config["params_annotate_snps"]["filter_nz_OTH"],
            filter_hom_ALT=config["params_annotate_snps"]["filter_hom_ALT"],
        threads: 1
        log:
            config["log_dir"] + "/annotate_snps_pseudobulk.log",
        script:
            "../scripts/annotate_snps_pseudobulk.py"
