##################################################
# SNP-informed adaptive binning + window→bin depth aggregation
# Bulk: combine_counts (with corrected window depth)
# Non-bulk: combine_counts_nonbulk + cnv_segmentation
##################################################


rule combine_counts:
    input:
        dp_corrected=lambda wc: config["pileup_dir"] + f"/{wc.assay_type}/window.dp.npz",
        window_df=lambda wc: config["pileup_dir"] + f"/{wc.assay_type}/window.tsv.gz",
        snp_info=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snps.tsv.gz",
        tot_mtx_snp=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/snp.Tallele.npz",
        a_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Aallele.npz",
        b_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Ballele.npz",
        sample_file=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
        gmap_file=lambda wc: (
            config["phase_dir"] + "/genetic_map.tsv.gz" if require_genetic_map else []
        ),
        region_bed=lambda wc: config["region_bed"],
        blacklist_bed=lambda wc: config.get("blacklist_bed") or [],
        genome_size=lambda wc: config["genome_size"],
        gtf_file=lambda wc: config["gtf_file"],
    output:
        bb_file=config["bb_dir"] + "/{assay_type}/bb.tsv.gz",
        tot_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Tallele.npz",
        a_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Aallele.npz",
        b_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Ballele.npz",
        baf_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.baf.npz",
        dp_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.depth.npz",
        rdr_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.rdr.npz",
        sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWES)",
    params:
        qc_dir=lambda wc: config["qc_dir"] + f"/{wc.assay_type}/combine_counts/",
        assay_type=lambda wc: wc.assay_type,
        nu=config["params_combine_counts"]["nu"],
        min_switchprob=config["params_combine_counts"]["min_switchprob"],
        switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
        min_snp_reads=config["params_combine_counts"]["min_snp_reads"],
        min_snp_per_block=config["params_combine_counts"]["min_snp_per_block"],
        skip_normal_normalization=config["params_combine_counts"].get("skip_normal_normalization", False),
        rdr_outlier_quantile=config["params_combine_counts"].get("rdr_outlier_quantile", 0.002),
        max_blocksize=config["params_combine_counts"].get("max_blocksize", 500000),
        phase_flip_test=config["params_combine_counts"].get("phase_flip_test", False),
        phase_flip_epsilon=config["params_combine_counts"].get("phase_flip_epsilon", 0.05),
        phase_flip_alpha=config["params_combine_counts"].get("phase_flip_alpha", 0.05),
        chromosomes=config["chromosomes"],
        run_id=_run_id,
    threads: 1
    log:
        config["log_dir"] + f"/combine_counts.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/combine_counts.py"""


rule combine_counts_nonbulk:
    input:
        snp_info=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snps.tsv.gz",
        tot_mtx_snp=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/snp.Tallele.npz",
        a_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Aallele.npz",
        b_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Ballele.npz",
        sample_file=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
        all_barcodes=lambda wc: branch(
            wc.assay_type in nonbulk_assays,
            then=config["allele_dir"] + f"/{wc.assay_type}/barcodes.tsv.gz",
            otherwise=[],
        ),
        gmap_file=lambda wc: (
            config["phase_dir"] + "/genetic_map.tsv.gz" if require_genetic_map else []
        ),
        region_bed=lambda wc: config["region_bed"],
        genome_size=lambda wc: config["genome_size"],
        gtf_file=lambda wc: config["gtf_file"],
    output:
        multi_snp_file=config["bb_dir"] + "/{assay_type}/multi_snp.tsv.gz",
        tot_mtx_multi=config["bb_dir"] + "/{assay_type}/multi_snp.Tallele.npz",
        a_mtx_multi=config["bb_dir"] + "/{assay_type}/multi_snp.Aallele.npz",
        b_mtx_multi=config["bb_dir"] + "/{assay_type}/multi_snp.Ballele.npz",
        bb_file=config["bb_dir"] + "/{assay_type}/bb.tsv.gz",
        tot_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Tallele.npz",
        a_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Aallele.npz",
        b_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Ballele.npz",
        baf_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.baf.npz",
        sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    params:
        qc_dir=lambda wc: config["qc_dir"] + f"/{wc.assay_type}/combine_counts/",
        assay_type=lambda wc: wc.assay_type,
        nu=config["params_combine_counts"]["nu"],
        min_switchprob=config["params_combine_counts"]["min_switchprob"],
        switchprob_ps=config["params_combine_counts"]["switchprob_ps"],
        nsnp_multi=config["params_combine_counts"]["nsnp_multi"],
        min_snp_reads=config["params_combine_counts"]["min_snp_reads"],
        min_snp_per_block=config["params_combine_counts"]["min_snp_per_block"],
        run_id=_run_id,
    threads: 1
    log:
        config["log_dir"] + f"/combine_counts_nonbulk.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/combine_counts_nonbulk.py"""


rule cnv_segmentation:
    input:
        snp_info=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snps.tsv.gz",
        tot_mtx_snp=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/snp.Tallele.npz",
        a_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Aallele.npz",
        b_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Ballele.npz",
        sample_file=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
        all_barcodes=config["allele_dir"] + "/{assay_type}/barcodes.tsv.gz",
        h5ad_file=config["bb_dir"] + "/{assay_type}/{assay_type}.h5ad",
        region_bed=lambda wc: config["region_bed"],
        genome_size=lambda wc: config["genome_size"],
        gtf_file=lambda wc: config["gtf_file"],
        bb_file=lambda wc: config["bb_file"],
    output:
        cnv_segments=config["bb_dir"] + "/{assay_type}/cnv_segments.tsv",
        x_count=config["bb_dir"] + "/{assay_type}/bb.Xcount.npz",
        tot_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Tallele.npz",
        a_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Aallele.npz",
        b_mtx_bb=config["bb_dir"] + "/{assay_type}/bb.Ballele.npz",
        barcodes_out=config["bb_dir"] + "/{assay_type}/barcodes.tsv.gz",
        sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    params:
        qc_dir=lambda wc: config["qc_dir"] + f"/{wc.assay_type}/cnv_segmentation/",
        sample_name=SAMPLE_ID,
        assay_type=lambda wc: wc.assay_type,
        feature_type=lambda wc: assay_type2feature_type[wc.assay_type],
        run_id=_run_id,
    threads: 1
    log:
        config["log_dir"] + f"/cnv_segmentation.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/cnv_segmentation.py"""
