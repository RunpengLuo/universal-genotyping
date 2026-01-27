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
elif config["phaser"] == "eagle":
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


rule concat_phased_snps:
    input:
        vcf_files=expand("phase/chr{chrname}.vcf.gz", chrname=config["chromosomes"]),
    output:
        phased_vcf="phase/phased.vcf.gz",
    threads: 1
    params:
        bcftools=config["bcftools"],
    shell:
        r"""
        printf "%s\n" {input.vcf_files} > "phase/phased_snps.lst"
        {params.bcftools} concat -f "phase/phased_snps.lst" -Oz -o "{output.phased_vcf}"
        tabix -f -p vcf "{output.phased_vcf}"
        rm "phase/phased_snps.lst"
        """

rule extract_phased_het_snps:
    input:
        phased_vcf="phase/phased.vcf.gz"
    output:
        phased_het_vcf="phase/phased_snps.vcf.gz",
    threads: 1
    params:
        bcftools=config["bcftools"],
    shell:
        r"""
        {params.bcftools} view "{input.phased_vcf}" -Oz -m2 -M2 \
            -i 'GT="0|1" || GT="1|0"' -o "{output.phased_het_vcf}"
        tabix -f -p vcf "{output.phased_het_vcf}"
        """
