##################################################
# Read depth computation and bias correction (bulk only)
#
# bulkWGS: mosdepth fixed-window depth + HMMcopy-style LOWESS correction
# bulkWES: CNVkit autobin/coverage/reference/fix pipeline
#
# Both produce window.dp.npz + window.tsv.gz consumed by combine_counts.
##################################################

_rdr_cfg = config.get("params_count_reads", {})


rule window_bed_to_3bed:
    """Extract headerless 3-column BED from window_bed for mosdepth --by."""
    input:
        window_bed=config.get("window_bed") or [],
    output:
        mosdepth_bed=temp(config["bb_dir"] + "/windows.bed.gz"),
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
        windows_bed=config["bb_dir"] + "/windows.bed.gz",
    output:
        mosdepth_file=config["bb_dir"]
        + "/{assay_type}/out_mosdepth/{rep_id}.regions.bed.gz",
    threads: config["threads"]["mosdepth"]
    wildcard_constraints:
        assay_type="bulkWGS",
    params:
        mosdepth=config["mosdepth"],
        out_prefix=config["bb_dir"] + "/{assay_type}/out_mosdepth/{rep_id}",
        read_quality=config["params_mosdepth"]["read_quality"],
        extra_params=config["params_mosdepth"].get("extra_params", ""),
    log:
        config["log_dir"] + f"/run_mosdepth/run_mosdepth.{{assay_type}}_{{rep_id}}.{_run_id}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        r"""
        {params.mosdepth} \
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
            config["bb_dir"] + f"/{wc.assay_type}/out_mosdepth/{rep_id}.regions.bed.gz"
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
        dp_corrected=config["bb_dir"] + "/{assay_type}/window.dp.npz",
        window_df=config["bb_dir"] + "/{assay_type}/window.tsv.gz",
        qc_dir=directory(config["qc_dir"] + "/{assay_type}/rd_correction/"),
    wildcard_constraints:
        assay_type="bulkWGS",
    params:
        sample_name=SAMPLE_ID,
        mosdepth_dir=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/out_mosdepth",
        chromosomes=config["chromosomes"],
        samplesize=_rdr_cfg.get("samplesize", 50000),
        routlier=_rdr_cfg.get("routlier", 0.01),
        doutlier=_rdr_cfg.get("doutlier", 0.001),
        min_mappability=_rdr_cfg.get("min_mappability", 0.9),
        gc_correct=_rdr_cfg.get("gc_correct", True),
        gc_correct_method=_rdr_cfg.get("gc_correct_method", "lowess"),
        rt_correct=_rdr_cfg.get("rt_correct", False),
        assay_type=lambda wc: wc.assay_type,
        run_id=_run_id,
    log:
        config["log_dir"] + f"/rd_correct.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/rd_correct.py"""


_cnvkit_cfg = config.get("params_cnvkit", {})
_cnvkit_exe = config.get("cnvkit", "cnvkit.py")
_cnvkit_threads = config.get("threads", {}).get("cnvkit", 4)


rule cnvkit_autobin:
    """Estimate bin sizes from a representative BAM; produce target/antitarget BEDs."""
    input:
        bam=lambda wc: branch(
            len(normal_bams) > 0, then=normal_bams[0], otherwise=tumor_bams[0]
        ),
        targets_bed=config.get("wes_targets_bed") or [],
        access_bed=config["region_bed"],
    output:
        target_bed=config["bb_dir"] + "/{assay_type}/cnvkit/targets.bed",
        antitarget_bed=config["bb_dir"] + "/{assay_type}/cnvkit/antitargets.bed",
    wildcard_constraints:
        assay_type="bulkWES",
    params:
        cnvkit=_cnvkit_exe,
        reference=config["reference"],
        chroms="|".join(f"chr{c}" for c in config["chromosomes"]),
        extra_flags=cli_flags_str(
            _cnvkit_cfg,
            ("method", "--method"),
            ("bp_per_bin", "--bp-per-bin"),
            ("target_max_size", "--target-max-size"),
            ("target_min_size", "--target-min-size"),
            ("antitarget_max_size", "--antitarget-max-size"),
            ("antitarget_min_size", "--antitarget-min-size"),
            ("annotate", "--annotate"),
            ("short_names", "--short-names", True),
        ),
    log:
        config["log_dir"] + f"/cnvkit_autobin.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        r"""
        {params.cnvkit} autobin {input.bam} \
            -t {input.targets_bed} \
            -g {input.access_bed} \
            -f {params.reference} \
            --target-output-bed {output.target_bed} \
            --antitarget-output-bed {output.antitarget_bed} \
            {params.extra_flags} > {log} 2>&1
        # Filter to standard chromosomes (remove alt/random/Un contigs)
        grep -E '^({params.chroms})\b' {output.target_bed} > {output.target_bed}.tmp \
            && mv {output.target_bed}.tmp {output.target_bed}
        grep -E '^({params.chroms})\b' {output.antitarget_bed} > {output.antitarget_bed}.tmp \
            && mv {output.antitarget_bed}.tmp {output.antitarget_bed}
        """


rule cnvkit_coverage:
    """Compute read coverage over target or antitarget regions per replicate."""
    input:
        bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
        interval=config["bb_dir"] + "/{assay_type}/cnvkit/{region_type}s.bed",
    output:
        cnn=config["bb_dir"] + "/{assay_type}/cnvkit/{rep_id}.{region_type}coverage.cnn",
    wildcard_constraints:
        assay_type="bulkWES",
        region_type="(target|antitarget)",
    params:
        cnvkit=_cnvkit_exe,
        extra_flags=cli_flags_str(_cnvkit_cfg, ("min_mapq", "-q")),
    threads: _cnvkit_threads
    log:
        config["log_dir"] + f"/cnvkit_coverage.{{assay_type}}_{{rep_id}}.{{region_type}}.{_run_id}.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        r"""
        {params.cnvkit} coverage {input.bam} {input.interval} \
            -o {output.cnn} \
            -p {threads} \
            {params.extra_flags} > {log} 2>&1
        """


rule cnvkit_reference:
    """Build a pooled reference from normal-sample coverage files."""
    input:
        normal_cnn=lambda wc: [
            config["bb_dir"] + f"/{wc.assay_type}/cnvkit/{rep_id}.{rt}coverage.cnn"
            for rep_id, st in zip(
                assay2rep_ids[wc.assay_type],
                assay2sample_types[wc.assay_type],
            )
            if st == "normal"
            for rt in ["target", "antitarget"]
        ],
    output:
        reference_cnn=config["bb_dir"] + "/{assay_type}/cnvkit/reference.cnn",
    wildcard_constraints:
        assay_type="bulkWES",
    params:
        cnvkit=_cnvkit_exe,
        reference=config["reference"],
        extra_flags=cli_flags_str(
            _cnvkit_cfg,
            ("no_gc", "--no-gc", True),
            ("no_edge", "--no-edge", True),
            ("no_rmask", "--no-rmask", True),
            ("male_reference", "--male-reference", True),
            ("cluster", "--cluster", True),
        ),
    log:
        config["log_dir"] + f"/cnvkit_reference.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        r"""
        {params.cnvkit} reference {input.normal_cnn} \
            -f {params.reference} \
            -o {output.reference_cnn} \
            {params.extra_flags} > {log} 2>&1
        """


rule cnvkit_fix:
    """Correct each replicate's coverage against the pooled reference."""
    input:
        target_cnn=config["bb_dir"] + "/{assay_type}/cnvkit/{rep_id}.targetcoverage.cnn",
        antitarget_cnn=config["bb_dir"]
        + "/{assay_type}/cnvkit/{rep_id}.antitargetcoverage.cnn",
        reference_cnn=config["bb_dir"] + "/{assay_type}/cnvkit/reference.cnn",
    output:
        cnr=config["bb_dir"] + "/{assay_type}/cnvkit/{rep_id}.cnr",
    wildcard_constraints:
        assay_type="bulkWES",
    params:
        cnvkit=_cnvkit_exe,
        extra_flags=cli_flags_str(
            _cnvkit_cfg,
            ("no_gc", "--no-gc", True),
            ("no_edge", "--no-edge", True),
            ("no_rmask", "--no-rmask", True),
        ),
    log:
        config["log_dir"] + f"/cnvkit_fix.{{assay_type}}_{{rep_id}}.{_run_id}.log",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        r"""
        {params.cnvkit} fix {input.target_cnn} {input.antitarget_cnn} {input.reference_cnn} \
            -o {output.cnr} \
            {params.extra_flags} > {log} 2>&1
        """


rule cnvkit_to_window_dp:
    """Convert .cnr files to window.dp.npz + window.tsv.gz for combine_counts."""
    input:
        cnr_files=lambda wc: [
            config["bb_dir"] + f"/{wc.assay_type}/cnvkit/{rep_id}.cnr"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        reference_cnn=config["bb_dir"] + "/{assay_type}/cnvkit/reference.cnn",
        sample_file=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
        region_bed=config["region_bed"],
        blacklist_bed=config.get("blacklist_bed") or [],
        genome_size=config["genome_size"],
    output:
        dp_corrected=config["bb_dir"] + "/{assay_type}/window.dp.npz",
        window_df=config["bb_dir"] + "/{assay_type}/window.tsv.gz",
        qc_dir=directory(config["qc_dir"] + "/{assay_type}/rd_correction/"),
    wildcard_constraints:
        assay_type="bulkWES",
    params:
        sample_name=SAMPLE_ID,
        cnvkit_dir=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/cnvkit",
        chromosomes=config["chromosomes"],
        run_id=_run_id,
    log:
        config["log_dir"] + f"/cnvkit_to_window_dp.{{assay_type}}.{_run_id}.log",
    conda:
        "../envs/base.yaml"
    script:
        """../scripts/cnvkit_to_window_dp.py"""
