##################################################
# Read depth computation and bias correction (bulk only)
#
# bulkWGS/bulkWES: mosdepth fixed-window depth + HMMcopy-style LOWESS correction
#
# Produces window.dp.npz + window.tsv.gz consumed by combine_counts.
##################################################

_rdr_cfg = config.get("params_count_reads", {})


rule window_bed_to_3bed:
    """Extract headerless 3-column BED from window_bed for mosdepth --by."""
    input:
        window_bed=config.get("window_bed") or [],
    output:
        mosdepth_bed=temp(config["pileup_dir"] + "/windows.bed.gz"),
    log:
        config["log_dir"] + f"/window_bed_to_3bed.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    run:
        import gzip, logging
        import pandas as pd

        logging.basicConfig(
            filename=str(log[0]),
            level=logging.INFO,
            format="%(asctime)s %(levelname)s %(message)s",
        )

        df = pd.read_table(
            str(input.window_bed), sep="\t", usecols=["#CHR", "START", "END"]
        )
        logging.info(f"Extracted {len(df)} windows from {input.window_bed}")
        with gzip.open(str(output.mosdepth_bed), "wt") as fh:
            df.to_csv(fh, sep="\t", header=False, index=False)


rule run_mosdepth:
    input:
        bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
        windows_bed=config["pileup_dir"] + "/windows.bed.gz",
    output:
        mosdepth_file=config["pileup_dir"]
        + "/{assay_type}/out_mosdepth/{rep_id}.regions.bed.gz",
    threads: config["threads"]["mosdepth"]
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWES)",
    params:
        out_prefix=config["pileup_dir"] + "/{assay_type}/out_mosdepth/{rep_id}",
        read_quality=config["params_mosdepth"]["read_quality"],
        extra_params=config["params_mosdepth"].get("extra_params", ""),
    log:
        config["log_dir"] + f"/run_mosdepth/run_mosdepth.{{assay_type}}_{{rep_id}}.{_run_id}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        r"""
        mosdepth \
            -t {threads} \
            -Q {params.read_quality} \
            --by {input.windows_bed} \
            {params.extra_params} \
            {params.out_prefix} {input.bam} > {log} 2>&1
        """


rule rd_correct:
    """Per-window LOWESS bias correction (GC/mappability/replication timing)."""
    input:
        mosdepth_files=lambda wc: [
            config["pileup_dir"] + f"/{wc.assay_type}/out_mosdepth/{rep_id}.regions.bed.gz"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        sample_file=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
        window_bed=config.get("window_bed") or [],
        genome_size=config["genome_size"],
        region_bed=config["region_bed"],
        blacklist_bed=config.get("blacklist_bed") or [],
        snp_info=config["allele_dir"] + "/{assay_type}/snps.tsv.gz",
        tot_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Tallele.npz",
        a_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Aallele.npz",
        b_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Ballele.npz",
    output:
        dp_corrected=config["pileup_dir"] + "/{assay_type}/window.dp.npz",
        window_df=config["pileup_dir"] + "/{assay_type}/window.tsv.gz",
        qc_dir=directory(config["qc_dir"] + "/{assay_type}/rd_correction/"),
    wildcard_constraints:
        assay_type="(bulkWGS|bulkWES)",
    params:
        sample_name=SAMPLE_ID,
        mosdepth_dir=lambda wc: config["pileup_dir"] + f"/{wc.assay_type}/out_mosdepth",
        chromosomes=config["chromosomes"],
        samplesize=_rdr_cfg.get("samplesize", 50000),
        routlier=_rdr_cfg.get("routlier", 0.01),
        doutlier=_rdr_cfg.get("doutlier", 0.001),
        min_mappability=_rdr_cfg.get("min_mappability", 0.9),
        gc_correct=_rdr_cfg.get("gc_correct", True),
        gc_correct_method=_rdr_cfg.get("gc_correct_method", "median"),
        rt_correct=_rdr_cfg.get("rt_correct", False),
        assay_type=lambda wc: wc.assay_type,
        run_id=_run_id,
    log:
        config["log_dir"] + f"/rd_correct.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/rd_correct.py"""


