##################################################
if config["phaser"] == "shapeit":
    rule phase_snps_shapeit:
        input:
            chrom_vcf_file=lambda wc: f"snps/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_gmap_file(wc.chrname),
        output:
            phased_file="phase/chr{chrname}.vcf.gz",
            bcf_file=temp("phase/chr{chrname}.bcf"),
            bcf_file_csi=temp("phase/chr{chrname}.bcf.csi"),
        threads: config["threads"]["phase"]
        params:
            chrom="chr{chrname}",
            shapeit=config["shapeit"],
            bcftools=config["bcftools"],
        log:
            "logs/phase_snps.chr{chrname}.log"
        shell:
            r"""
            {params.shapeit} \
                --input "{input.chrom_vcf_file}" \
                --map "{input.gmap_file}" \
                --reference "{input.phasing_panel_file}" \
                --region "{params.chrom}" \
                --thread "{threads}" \
                --output "{output.bcf_file}"
            
            bcftools view -Ov "{output.bcf_file}" | bgzip > "{output.phased_file}"
            tabix -f -p vcf "{output.phased_file}"
            """

##################################################
if config["phaser"] == "eagle":
    rule phase_snps_eagle:
        input:
            chrom_vcf_file=lambda wc: f"snps/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_gmap_file(wc.chrname),
        output:
            phased_file="phase/chr{chrname}.vcf.gz",
        threads: config["threads"]["phase"]
        params:
            chrom="chr{chrname}",
            eagle=config["eagle"],
            bcftools=config["bcftools"],
        log:
            "logs/phase_snps.chr{chrname}.log"
        shell:
            r"""
            {params.eagle} \
                --vcfTarget "{input.chrom_vcf_file}" \
                --geneticMapFile "{input.gmap_file}" \
                --vcfRef "{input.phasing_panel_file}" \
                --vcfOutFormat z \
                --numThreads "{threads}" \
                --outPrefix "phase/{params.chrom}"
            tabix -f -p vcf "{output.phased_file}"
            """

##################################################
if config["phaser"] == "longphase":
    rule phase_snps_longphase:
        input:
            chrom_vcf_file=lambda wc: f"snps/chr{wc.chrname}.vcf.gz",
            bam_file=lambda wc: bulk_nbams[0] if has_normal else bulk_tbams[0],
            ref_fa=lambda wc: config["reference"]
        output:
            phased_file="phase/chr{chrname}.vcf.gz",
        params:
            chrom="chr{chrname}",
            longphase=config["longphase"],
            min_mapq=config["params_longphase"]["min_mapq"],
            extra_params=config["params_longphase"].get("extra_params", ""),
            bcftools=config["bcftools"],
        threads: config["threads"]["phase"]
        log: 
            "logs/phase_snps.chr{chrname}.log"
        shell:
            r"""
            {params.longphase} phase \
                --bam-file={input.bam_file} \
                --reference={input.ref_fa} \
                --snp-file={input.chrom_vcf_file} \
                --mappingQuality={params.min_mapq} \
                --out-prefix="phase/{params.chrom}" \
                --threads={threads} \
                {params.extra_params}
            tabix -f -p vcf "{output.phased_file}"
            """

rule concat_and_extract_phased_het_snps:
    input:
        vcf_files=expand("phase/chr{chrname}.vcf.gz", chrname=config["chromosomes"]),
    output:
        phased_vcf="phase/phased_snps.vcf.gz",
        lst_file=temp("phase/phased_snps.lst")
    threads: 1
    params:
        bcftools=config["bcftools"],
    shell:
        r"""
        printf "%s\n" {input.vcf_files} > "{output.lst_file}"
        {params.bcftools} concat -f "{output.lst_file}" -Ou \
        | {params.bcftools} view -Oz -m2 -M2 -i 'GT="0|1" || GT="1|0"' \
            -o "{output.phased_vcf}"
        tabix -f -p vcf "{output.phased_vcf}"
        """

rule parse_genetic_map:
    input:
        gmap_files=lambda wc: [get_gmap_file(c) for c in config["chromosomes"]],
    output:
        gmap_tsv="phase/genetic_map.tsv.gz",
    log:
        "logs/parse_genetic_map.log",
    params:
        chrnames=config["chromosomes"],
        phaser=config["phaser"],
    threads: 1
    script:
        "../scripts/parse_genetic_map.py"
        