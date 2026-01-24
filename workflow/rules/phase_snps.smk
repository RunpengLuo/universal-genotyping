if config["phaser"] == "shapeit":
    rule phase_snps_shapeit:
        input:
            chrom_vcf_file=lambda wc: f"snps/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_gmap_file(wc.chrname),
        output:
            done_file="phase/.shapeit.chr{chrname}.done",
            phased_file="phase/phased.chr{chrname}.vcf.gz",
        threads: config["threads"]["phase"]
        params:
            chrom="chr{chrname}",
            shapeit=config["shapeit"],
            bcftools=config["bcftools"],
        log:
            "logs/phase_snps.chr{chrname}.log"
        shell:
            r"""
            set -euo pipefail
            {params.shapeit} \
                --input "{input.chrom_vcf_file}" \
                --map "{input.gmap_file}" \
                --reference "{input.phasing_panel_file}" \
                --region "{params.chrom}" \
                --thread "{threads}" \
                --output "phase/phased.{params.chrom}.vcf" >> "{log}" 2>&1
            
            bgzip -f "phase/phased.{params.chrom}.vcf"
            tabix -f -p vcf "{output.phased_file}"
            touch "{output.done_file}"
            """
elif config["phaser"] == "eagle":
    rule phase_snps_eagle:
        input:
            chrom_vcf_file=lambda wc: f"snps/chr{wc.chrname}.vcf.gz",
            phasing_panel_file=lambda wc: get_phasing_panel(wc.chrname),
            gmap_file=lambda wc: get_gmap_file(wc.chrname),
        output:
            done_file="phase/.eagle.chr{chrname}.done",
            phased_file="phase/phased.chr{chrname}.vcf.gz",
        threads: config["threads"]["phase"]
        params:
            chrom="chr{chrname}",
            eagle=config["eagle"],
            bcftools=config["bcftools"],
        log:
            "logs/phase_snps.chr{chrname}.log"
        shell:
            r"""
            set -euo pipefail
            {params.eagle} \
                --vcfTarget "{input.chrom_vcf_file}" \
                --geneticMapFile "{input.gmap_file}" \
                --vcfRef "{input.phasing_panel_file}" \
                --vcfOutFormat z \
                --numThreads "{threads}" \
                --outPrefix "phase/phased.{params.chrom}" >> "{log}" 2>&1
            tabix -f -p vcf "{output.phased_file}"
            touch "{output.done_file}"      
            """


rule concat_phased_snps:
    input:
        vcf_files=expand("phase/phased.chr{chrname}.vcf.gz", chrname=config["chromosomes"]),
    output:
        phased_vcf="phase/phased.vcf.gz",
        lst_file=temp("phase/phased_snps.lst"),
    threads: 1
    params:
        bcftools=config["bcftools"],
    shell:
        r"""
        set -euo pipefail
        printf "%s\n" {input.vcf_files} > "{output.lst_file}"
        {params.bcftools} concat -f "{output.lst_file}" -Ou \
        | {params.bcftools} view -g het -Oz -o "{output.phased_vcf}"
        tabix -f -p vcf "{output.phased_vcf}"
        """

rule extract_phased_het_snps:
    input:
        phased_vcf="phase/phased.vcf.gz"
    output:
        phased_het_vcf="phase/phased_snps.vcf.gz",
        phased_het_vcfgz="phase/phased_snps.vcf.gz.tbi"
    threads: 1
    params:
        bcftools=config["bcftools"],
    shell:
        r"""
        set -euo pipefail
        {params.bcftools} view -Ou {input.phased_vcf} \
            -v snps -m2 -M2 -i 'GT="0|1" || GT="1|0"' \
        | {params.bcftools} view -Oz -o "{output.phased_het_vcf}"
        tabix -f -p vcf "{output.phased_het_vcf}"
        """
