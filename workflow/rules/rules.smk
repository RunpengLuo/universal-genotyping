from snakemake.io import directory

##################################################
rule genotype_snps:
    """
    Genotype SNPs from pseudobulk BAMs via cellsnp-lite mode 1b.
    """
    input:
        bams=bams,
        snp_panel=config["snp_panel"],
    output:
        raw_snp_file="genotype/cellSNP.base.vcf.gz",
    threads: config["threads"]["genotype"]
    params:
        cellsnp_lite=config["cellsnp_lite"],
        UMItag=config["params_cellsnp_lite"]["UMItag"],
        minMAF=config["params_cellsnp_lite"]["minMAF"],
        minCOUNT=config["params_cellsnp_lite"]["minCOUNT"],
    log:
        "logs/genotype_snps.log",
    shell:
        r"""
        set -euo pipefail
        echo "genotype SNPs on pseudobulk sample"
        printf "%s\n" {input.bams} > genotype/bams.lst
        {params.cellsnp_lite} \
            -S genotype/bams.lst \
            -O genotype \
            -R {input.snp_panel} \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag {params.UMItag} \
            --cellTAG None \
            --gzip > "{log}" 2>&1
        """

##################################################
rule annotate_snps:
    """
    Annotate annotate SNP GT information, 
    pre-filter SNP non-hom/het SNPs before phasing.
    """
    input:
        raw_snp_file="genotype/cellSNP.base.vcf.gz",
    output:
        expand("genotype/chr{chrname}.vcf", chrname=config["chromosomes"]),
    threads: 1
    log:
        "logs/annotate_snps.log"
    script:
        "../scripts/annotate_snps.py"

rule bgzip_index_snps:
    input: 
        vcf="genotype/chr{chrname}.vcf"
    output: 
        vcfgz="genotype/chr{chrname}.vcf.gz",
        tbi="genotype/chr{chrname}.vcf.gz.tbi"
    threads: 1
    shell: 
        """
        set -euo pipefail
        bgzip -f {input.vcf}
        tabix -f -p vcf {output.vcfgz}
        """

##################################################
rule phase_snps_per_chrom:
    """
    Phase genotyped SNPs via population-based phasing methods
    """
    input:
        chrom_vcf_file=lambda wc: f"genotype/chr{wc.chrname}.vcf.gz",
        phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
        gmap_file=lambda wc: get_gmap_file(wc.chrname),
    output:
        "phase/phased.chr{chrname}.vcf.gz"
    threads: config["nthreads"]["phase"]
    params:
        chrom="chr{chrname}",
        phaser=config["phaser"],
        eagle=config["eagle"],
        shapeit=config["shapeit"],
        bcftools=config["bcftools"],
    log:
        "logs/phase_snps.chr{chrname}.log"
    shell:
        r"""
        set -euo pipefail
        if [[ "{params.phaser}" == "shapeit" ]]; then
            {params.shapeit} \
                --input "{input.chrom_vcf_file}" \
                --map "{input.gmap_file}" \
                --reference "{input.phasing_panel_file}" \
                --region "{params.chrom}" \
                --thread "{threads}" \
                --output "phase/phased.{params.chrom}.vcf" > "{log}" 2>&1
            {params.bcftools} view -Oz \
                -o "{output}" "phase/phased.{params.chrom}.vcf"
            rm "phase/phased.{params.chrom}.vcf"
        elif [[ "{params.phaser}" == "eagle" ]]; then
            {params.eagle} \
                --vcfTarget "{input.chrom_vcf_file}" \
                --geneticMapFile "{input.gmap_file}" \
                --vcfRef "{input.phasing_panel_file}" \
                --vcfOutFormat z \
                --numThreads "{threads}" \
                --outPrefix "phase/phased.{params.chrom}" > "{log}" 2>&1
        else
            echo "undefined phase method: {params.phaser}" >&2
            exit 1
        fi
        tabix -f -p vcf {output}
        bcftools index -f "{output}"
        """

rule concat_phased_snps:
    input:
        expand("phase/phased.chr{chrname}.vcf.gz", chrname=config["chromosomes"])
    output:
        "phase/phased.vcf.gz"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        printf "%s\n" {input} > phase/concat.lst
        bcftools concat -f phase/concat.lst -Ou \
          | bcftools view -g het -Oz -o "{output}"
        bcftools index -f "{output}"
        """


##################################################
rule pileup_snps:
    """
    Pile-up cell/spot by feature allele count matrices via cellsnp-lite mode1a
    """
    input:
        barcode=lambda wc: get_files[(wc.mod, wc.rep_id)][0],
        bam=lambda wc: get_files[(wc.mod, wc.rep_id)][1],
        ranger=lambda wc: get_files[(wc.mod, wc.rep_id)][2],
        snp_file=lambda wc: ("phase/phased.vcf.gz" if require_genotyping else ref_snp_file),
    output:
        out_dir=directory("pileup/{mod}_{rep_id}"),
        cellsnp_file="pileup/{mod}_{rep_id}/cellSNP.base.vcf.gz",
    threads: config["threads"]["pileup"]
    params:
        cellsnp_lite=config["cellsnp_lite"],
        UMItag=lambda wc: (config["params_cellsnp_lite"]["UMItag"] if wc.mod != "ATAC" else "None"),
        cellTAG=config["params_cellsnp_lite"]["cellTAG"],
        minMAF=config["params_cellsnp_lite"]["minMAF"],
        minCOUNT=config["params_cellsnp_lite"]["minCOUNT"],
    log:
        "logs/pileup_snps.{mod}_{rep_id}.log"
    shell:
        r"""
        set -euo pipefail
        echo "pileup SNPs, out_dir={output.out_dir}"
        mkdir -p "{output.out_dir}"

        {params.cellsnp_lite} \
            -b "{input.barcode}" \
            -s "{input.bam}" \
            -O "{output.out_dir}" \
            -R "{input.snp_file}" \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT} \
            --UMItag {params.UMItag} \
            --cellTAG {params.cellTAG} \
            --gzip > "{log}" 2>&1
        """

##################################################
rule postprocess_matrix:
    input:
        vcfs=expand("pileup/{mod}_{rep_id}/cellSNP.base.vcf.gz", mod=mods, rep_id=rep_ids),
        snp_file=("phase/phased.vcf.gz" if require_genotyping else ref_snp_file)
    output:
        bc_mod=expand("{mod}/barcodes.txt", mod=mods),
        A_mod=expand("{mod}/cell_snp_Aallele.npz", mod=mods),
        B_mod=expand("{mod}/cell_snp_Ballele.npz", mod=mods),
    threads: config["threads"]["postprocess"]
    params:
        mods=mods,
        rep_ids=rep_ids,
    log:
        "logs/postprocess_matrix.log"
    script:
        "../scripts/postprocess_matrix.py"
