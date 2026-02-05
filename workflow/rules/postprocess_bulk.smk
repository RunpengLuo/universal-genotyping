##################################################
if workflow_mode == "bulk":

    rule postprocess_bulk:
        input:
            vcfs=lambda wc: expand(
                "pileup/bulkDNA_{rep_id}/cellSNP.base.vcf.gz",
                rep_id=mod2reps["bulkDNA"],
            ),
            tot_mats=lambda wc: expand(
                "pileup/bulkDNA_{rep_id}/cellSNP.tag.DP.mtx",
                rep_id=mod2reps["bulkDNA"],
            ),
            ad_mats=lambda wc: expand(
                "pileup/bulkDNA_{rep_id}/cellSNP.tag.AD.mtx",
                rep_id=mod2reps["bulkDNA"],
            ),
            snp_file=lambda wc: "phase/phased_snps.vcf.gz",
            genome_size=lambda wc: config["genome_size"],
            region_bed=lambda wc: config["region_bed"],
            gmap_file=lambda wc: (
                "phase/genetic_map.tsv.gz" if require_genetic_map else []
            ),
            block_bed=lambda wc: (
                [] if config.get("block_bed") is None else config["block_bed"]
            ),
        output:
            sample_file="allele/bulkDNA/sample_ids.tsv",
            snp_file="allele/bulkDNA/snps.tsv.gz",
            tot_mtx_snp="allele/bulkDNA/snp.Tallele.npz",
            a_mtx_snp="allele/bulkDNA/snp.Aallele.npz",
            b_mtx_snp="allele/bulkDNA/snp.Ballele.npz",
            bb_file="allele/bulkDNA/bb.tsv.gz",
            bed_file="allele/bulkDNA/bb.bed.gz",
            tot_mtx_bb="allele/bulkDNA/bb.Tallele.npz",
            a_mtx_bb="allele/bulkDNA/bb.Aallele.npz",
            b_mtx_bb="allele/bulkDNA/bb.Ballele.npz",
            qc_dir=directory("allele/bulkDNA/qc"),
        params:
            sample_name=SAMPLE_ID,
            modality="bulkDNA",
            data_types=["bulkDNA"],
            rep_ids=mod2reps.get("bulkDNA", []),
            min_depth=config["params_postprocess"]["min_depth"],
            gamma=config["params_postprocess"]["gamma"],
            nu=config["params_postprocess"]["nu"],
            min_switchprob=config["params_postprocess"]["min_switchprob"],
            max_switchprob=config["params_postprocess"]["max_switchprob"],
            switchprob_ps=config["params_postprocess"]["switchprob_ps"],
            binom_test=config["params_postprocess"]["binom_test"],
            binom_alpha=config["params_postprocess"]["binom_alpha"],
            binom_margin=config["params_postprocess"]["binom_margin"],
            min_snp_reads=config["params_postprocess"]["min_snp_reads"],
            min_snp_per_block=config["params_postprocess"]["min_snp_per_block"],
        log:
            "logs/postprocess.bulkDNA.log",
        script:
            """../scripts/postprocess_bulk.py"""

    rule run_mosdepth_bulk:
        input:
            bam=lambda wc: get_data[("bulkDNA", wc.rep_id)][1],
            bed_file="allele/bulkDNA/bb.bed.gz",
        output:
            mosdepth_file="allele/bulkDNA/out_mosdepth/{rep_id}.regions.bed.gz",
        threads: config["threads"]["mosdepth"]
        params:
            mosdepth=config["mosdepth"],
            sample_name=SAMPLE_ID,
            modality="bulkDNA",
            out_prefix="allele/bulkDNA/out_mosdepth/{rep_id}",
            read_quality=config["params_mosdepth"]["read_quality"],
            extra_params=config["params_mosdepth"].get("extra_params", ""),
        log:
            "logs/mosdepth.bulkDNA_{rep_id}.log",
        shell:
            r"""
            {params.mosdepth} \
                -t {threads} \
                -Q {params.read_quality} \
                --by "{input.bed_file}" \
                {params.extra_params} \
                "{params.out_prefix}" \
                "{input.bam}" > {log} 2>&1
            """

    rule compute_rdr_bulk:
        input:
            sample_file="allele/bulkDNA/sample_ids.tsv",
            bb_file="allele/bulkDNA/bb.tsv.gz",
            mosdepth_files=expand(
                "allele/bulkDNA/out_mosdepth/{rep_id}.regions.bed.gz",
                rep_id=mod2reps.get("bulkDNA", []),
            ),
            reference=config["reference"],
            genome_size=config["genome_size"],
            mappability_file=branch(
                config.get("mappability_file") is None,
                then=[],
                otherwise=config["mappability_file"],
            ),
        output:
            dp_mtx_bb="allele/bulkDNA/bb.depth.npz",
            rdr_mtx_bb="allele/bulkDNA/bb.rdr.npz",
        params:
            sample_name=SAMPLE_ID,
            qc_dir="allele/bulkDNA/qc",
            mosdepth_dir="allele/bulkDNA/out_mosdepth",
            gc_correct=config["params_compute_rdr"]["gc_correct"],
        log:
            "logs/compute_rdrs.bulkDNA.log",
        script:
            """../scripts/compute_rdr_bulk.py"""
