##################################################
if workflow_mode == "bulk":

    # TODO genotyping tumor samples
    rule genotype_snps_bulk:
        """
        Genotype Bi-allelic HET/ALT-HOM SNPs
        from bulk-DNA normal sample via bcftools
        """
        input:
            bams=lambda wc: branch(
                has_normal, then=bulk_nbams[0], otherwise=bulk_tbams[0]
            ),
            snp_panel=config["snp_panel"],
            reference=config["reference"],
        output:
            snp_file="snps/chr{chrname}.vcf.gz",
            tmp_pos=temp("tmp/target.{chrname}.pos.gz"),
            tmp_pos_tbi=temp("tmp/target.{chrname}.pos.gz.tbi"),
        log:
            "logs/genotype_snps.chr{chrname}.log",
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
                -Oz -o {output.snp_file} > {log} 2>&1

            tabix -p vcf {output.snp_file}
            """


##################################################
if workflow_mode == "single_cell":
    rule genotype_snps_pseudobulk:
        input:
            bams=gex_tbams+atac_tbams,
            snp_panel=config["snp_panel"],
        output:
            snp_file="pileup/pseudobulk/cellSNP.base.vcf.gz",
        log:
            "logs/genotype_snps.pseudobulk.log"
        threads:
            config["threads"]["genotype"]
        params:
            bcftools=config["bcftools"],
            cellsnp_lite=config["cellsnp_lite"],
            out_dir=lambda wc: f"pileup/pseudobulk",
            minMAF=0,
            minCOUNT=1,
        shell:
            r"""
            printf "%s\n" {input.bams} > {params.out_dir}/bams.lst
            {params.cellsnp_lite} \
                -S {params.out_dir}/bams.lst \
                -O "{params.out_dir}" \
                -R "{input.snp_panel}" \
                -p {threads} \
                --minMAF {params.minMAF} \
                --minCOUNT {params.minCOUNT} \
                --UMItag None \
                --cellTAG None \
                --gzip > {log} 2>&1
            nsnps_cellsnp=$({params.bcftools} view -H "{output.snp_file}" | wc -l)
            echo "[QC] {output} record #SNPs: ${{nsnps_cellsnp}}"
            rm {params.out_dir}/bams.lst
            """

    rule genotype_snps_single_cell:
        input:
            barcode=lambda wc: get_data[(wc.data_type, wc.rep_id)][0],
            bam=lambda wc: get_data[(wc.data_type, wc.rep_id)][1],
            ranger=lambda wc: get_data[(wc.data_type, wc.rep_id)][2],
            snp_file="pileup/pseudobulk/cellSNP.base.vcf.gz",
        output:
            cellsnp_file="pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz",
            sample_file="pileup/{data_type}_{rep_id}/cellSNP.samples.tsv",
            tot_mat="pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx",
            ad_mat="pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx",
        wildcard_constraints:
            data_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        threads: config["threads"]["genotype"]
        params:
            cellsnp_lite=config["cellsnp_lite"],
            refseq=config["reference"],
            out_dir=lambda wc: f"pileup/{wc.data_type}_{wc.rep_id}",
            UMItag=lambda wc: branch(
                wc.data_type == "scATAC",
                then="None",
                otherwise=config["params_cellsnp_lite"]["UMItag"],
            ),
            cellTAG=config["params_cellsnp_lite"]["cellTAG"],
            minMAF=config["params_cellsnp_lite"]["minMAF"],
            minCOUNT=config["params_cellsnp_lite"]["minCOUNT_genotype"],
        log:
            "logs/genotype_snps.{data_type}_{rep_id}.log",
        shell:
            r"""
            {params.cellsnp_lite} \
                -b "{input.barcode}" \
                -s "{input.bam}" \
                -O "{params.out_dir}" \
                -R "{input.snp_file}" \
                -p {threads} \
                --minMAF {params.minMAF} \
                --minCOUNT {params.minCOUNT} \
                --UMItag {params.UMItag} \
                --cellTAG {params.cellTAG} \
                --gzip > {log} 2>&1
            """

    rule annotate_snps_single_cell:
        input:
            raw_snp_files=[
                f"pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz"
                for (data_type, rep_id) in get_data.keys()
            ],
            sample_files=[
                f"pileup/{data_type}_{rep_id}/cellSNP.samples.tsv"
                for (data_type, rep_id) in get_data.keys()
            ],
            dp_files=[
                f"pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx"
                for (data_type, rep_id) in get_data.keys()
            ],
            ad_files=[
                f"pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx"
                for (data_type, rep_id) in get_data.keys()
            ],
        output:
            snp_files=expand("snps/chr{chrname}.vcf.gz", chrname=config["chromosomes"]),
            snp_files_tbi=expand(
                "snps/chr{chrname}.vcf.gz.tbi", chrname=config["chromosomes"]
            ),
        params:
            data_types=[data_type for (data_type, _) in get_data.keys()],
            rep_ids=[rep_id for (_, rep_id) in get_data.keys()],
            filter_nz_OTH=config["params_annotate_snps"]["filter_nz_OTH"],
            filter_hom_ALT=config["params_annotate_snps"]["filter_hom_ALT"],
            min_het_reads=config["params_annotate_snps"]["min_het_reads"],
            min_hom_dp=config["params_annotate_snps"]["min_hom_dp"],
            min_vaf_thres=config["params_annotate_snps"]["min_vaf_thres"],
        threads: 1
        log:
            "logs/annotate_snps.log",
        script:
            "../scripts/annotate_snps.py"
