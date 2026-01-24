if genotype_mode == "bulk_normal":
    rule genotype_snps_bulk_normal:
        """
        Genotype Bi-allelic HET/ALT-HOM SNPs 
        from bulk-DNA normal sample via bcftools
        """
        input:
            bams=lambda wc: bulk_nbams[0] if len(bulk_nbams) > 0 else [],
            snp_panel=config["snp_panel"],
            reference=config["reference"]
        output:
            snp_file="snps/chr{chrname}.vcf.gz",
            tmp_pos=temp("tmp/target.{chrname}.pos.gz"),
            tmp_pos_tbi=temp("tmp/target.{chrname}.pos.gz.tbi")
        log: "logs/genotype_snps.{chrname}.log",
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
            set -euo pipefail

            {params.bcftools} query -f '%CHROM\t%POS\n' -r {params.chrom} \
                {input.snp_panel} | bgzip -c > {output.tmp_pos}
            tabix -s1 -b2 -e2 {output.tmp_pos}

            {params.bcftools} mpileup {input.bams} \
                -f "{input.reference}" \
                -Ou \
                -a INFO/AD,AD,DP \
                --skip-indels \
                -q {params.min_mapq} \
                -Q {params.min_baseq} \
                -d {params.max_depth} \
                -T {output.tmp_pos} \
            | {params.bcftools} call -m -Ou \
            | {params.bcftools} view -v snps -m2 -M2 \
                -i 'GT="alt" && FMT/DP>={params.min_dp}' \
                -Oz -o {output.snp_file} >> "{log}" 2>&1

            tabix -p vcf {output.snp_file}
            """

elif genotype_mode == "bulk_tumor":
    # genotype SNPs from bulk tumor samples
    # TODO
    pass 
elif genotype_mode == "pseudobulk_tumor":
    rule genotype_snps_pseudobulk:
        input:
            bams=lambda wc: {"sc_GEX": gex_tbams, "sc_ATAC": atac_tbams}[wc.data_type],
            snp_panel=config["snp_panel"],
        output:
            cellsnp_dir=directory("snps/{data_type}"),
            snp_file="snps/{data_type}/cellSNP.base.vcf.gz",
            lst_file=temp("snps/{data_type}_bams.lst"),
        log: "logs/genotype_snps.{data_type}.log",
        threads: config["threads"]["genotype"],
        params:
            bcftools=config["bcftools"],
            cellsnp_lite=config["cellsnp_lite"],
            chroms=",".join(map(str, config["chromosomes"])),
            refseq=config["reference"],
            UMItag=lambda wc: ("None" if wc.data_type == "sc_ATAC" else config["params_cellsnp_lite"]["UMItag"]),
            minMAF=config["params_cellsnp_lite"]["minMAF"],
            minCOUNT=config["params_cellsnp_lite"]["minCOUNT_genotype"],
        shell:
            r"""
            set -euo pipefail
            echo "genotype SNPs on pseudobulk {data_type} sample"
            printf "%s\n" {input.bams} > {output.lst_file}
            mkdir -p {output.cellsnp_dir}

            nsnps_panel=$({params.bcftools} view -H {input.snp_panel} | wc -l)
            echo "[QC] {input.snp_panel} record #SNPs: ${nsnps_panel}" >> "{log}"
            {params.cellsnp_lite} \
                -S {output.lst_file} \
                -O {output.cellsnp_dir} \
                -R {input.snp_panel} \
                --chrom "{params.chroms}" \
                --refseq {params.refseq} \
                -p {threads} \
                --minMAF {params.minMAF} \
                --minCOUNT {params.minCOUNT} \
                --UMItag {params.UMItag} \
                --cellTAG None \
                --gzip >> "{log}" 2>&1
            
            nsnps_cellsnp=$({params.bcftools} view -H "{output.snp_file}" | wc -l)
            echo "[QC] {output} record #SNPs: ${nsnps_cellsnp}" >> "{log}"
            """

    rule annotate_snps_pseudobulk:
        input:
            raw_snp_files=expand("snps/{data_type}/cellSNP.base.vcf.gz", data_type=genotype_dtypes)
        output:
            snp_files=expand("snps/chr{chrname}.vcf.gz", chrname=config["chromosomes"]),
            snp_files_tbi=expand("snps/chr{chrname}.vcf.gz.tbi", chrname=config["chromosomes"]),
        params:
            data_types=genotype_dtypes,
            filter_nz_OTH=config["params_annotate_snps"]["filter_nz_OTH"]
        threads: 1
        log:
            "logs/annotate_snps.log"
        script:
            "../scripts/annotate_snps.py"

