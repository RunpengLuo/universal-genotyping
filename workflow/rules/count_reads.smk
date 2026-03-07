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
        window_bed=config["window_bed"],
    output:
        mosdepth_bed=temp(config["bb_dir"] + "/windows.bed.gz"),
    log:
        config["log_dir"] + "/window_bed_to_3bed.log",
    conda:
        "../envs/base.yaml"
    run:
        import gzip, logging
        import pandas as pd
        logging.basicConfig(filename=str(log[0]), level=logging.INFO,
                            format="%(asctime)s %(levelname)s %(message)s")

        df = pd.read_table(str(input.window_bed), sep="\t", usecols=["#CHR", "START", "END"])
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
        config["log_dir"]
        + "/run_mosdepth/run_mosdepth.{assay_type}_{rep_id}.log",
    conda:
        "../envs/tools.yaml"
    run:
        import subprocess

        cmd = [
            str(params.mosdepth),
            "-t", str(threads),
            "-Q", str(params.read_quality),
            "--by", str(input.windows_bed),
        ]
        if params.extra_params:
            cmd.extend(params.extra_params.split())
        cmd.extend([str(params.out_prefix), str(input.bam)])

        with open(str(log[0]), "w") as fh:
            subprocess.check_call(cmd, stdout=fh, stderr=subprocess.STDOUT)


rule rd_correct:
    """Per-window LOWESS bias correction (GC/mappability/replication timing)."""
    input:
        mosdepth_files=lambda wc: [
            config["bb_dir"]
            + f"/{wc.assay_type}/out_mosdepth/{rep_id}.regions.bed.gz"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        sample_file=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/sample_ids.tsv",
        window_bed=config["window_bed"],
        genome_size=config["genome_size"],
        region_bed=config["region_bed"],
        blacklist_bed=config.get("blacklist_bed", []),
    output:
        dp_corrected=config["bb_dir"] + "/{assay_type}/window.dp.npz",
        window_df=config["bb_dir"] + "/{assay_type}/window.tsv.gz",
        qc_dir=directory(config["qc_dir"] + "/{assay_type}/rd_correction/"),
    wildcard_constraints:
        assay_type="bulkWGS",
    params:
        sample_name=SAMPLE_ID,
        mosdepth_dir=lambda wc: config["bb_dir"]
        + f"/{wc.assay_type}/out_mosdepth",
        chromosomes=config["chromosomes"],
        samplesize=_rdr_cfg.get("samplesize", 50000),
        routlier=_rdr_cfg.get("routlier", 0.01),
        doutlier=_rdr_cfg.get("doutlier", 0.001),
        min_mappability=_rdr_cfg.get("min_mappability", 0.9),
        gc_correct=_rdr_cfg.get("gc_correct", True),
        rt_correct=_rdr_cfg.get("rt_correct", False),
        assay_type=lambda wc: wc.assay_type,
    log:
        config["log_dir"] + "/rd_correct.{assay_type}.log",
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
        targets_bed=config.get("wes_targets_bed", ""),
        access_bed=config["region_bed"],
    output:
        target_bed=config["bb_dir"] + "/{assay_type}/cnvkit/targets.bed",
        antitarget_bed=config["bb_dir"] + "/{assay_type}/cnvkit/antitargets.bed",
    wildcard_constraints:
        assay_type="bulkWES",
    params:
        cnvkit=_cnvkit_exe,
        reference=config["reference"],
        cfg=_cnvkit_cfg,
    log:
        config["log_dir"] + "/cnvkit_autobin.{assay_type}.log",
    conda:
        "../envs/cnvkit.yaml"
    run:
        import subprocess

        cmd = [
            str(params.cnvkit), "autobin",
            str(input.bam),
            "-t", str(input.targets_bed),
            "-g", str(input.access_bed),
            "-f", str(params.reference),
            "--target-output-bed", str(output.target_bed),
            "--antitarget-output-bed", str(output.antitarget_bed),
        ]
        cfg = params.cfg
        cmd += cli_flag(cfg, "method", "--method")
        cmd += cli_flag(cfg, "bp_per_bin", "--bp-per-bin")
        cmd += cli_flag(cfg, "target_max_size", "--target-max-size")
        cmd += cli_flag(cfg, "target_min_size", "--target-min-size")
        cmd += cli_flag(cfg, "antitarget_max_size", "--antitarget-max-size")
        cmd += cli_flag(cfg, "antitarget_min_size", "--antitarget-min-size")
        cmd += cli_flag(cfg, "annotate", "--annotate")
        cmd += cli_flag(cfg, "short_names", "--short-names", is_bool=True)

        with open(str(log[0]), "w") as fh:
            subprocess.check_call(cmd, stdout=fh, stderr=subprocess.STDOUT)


rule cnvkit_coverage:
    """Compute read coverage over target or antitarget regions per replicate."""
    input:
        bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
        interval=config["bb_dir"] + "/{assay_type}/cnvkit/{region_type}s.bed",
    output:
        cnn=config["bb_dir"]
        + "/{assay_type}/cnvkit/{rep_id}.{region_type}coverage.cnn",
    wildcard_constraints:
        assay_type="bulkWES",
        region_type="(target|antitarget)",
    params:
        cnvkit=_cnvkit_exe,
        reference=config["reference"],
        cfg=_cnvkit_cfg,
    threads: _cnvkit_threads
    log:
        config["log_dir"]
        + "/cnvkit_coverage.{assay_type}_{rep_id}.{region_type}.log",
    conda:
        "../envs/cnvkit.yaml"
    run:
        import subprocess

        cmd = [
            str(params.cnvkit), "coverage",
            str(input.bam),
            str(input.interval),
            "-o", str(output.cnn),
            "-p", str(threads),
        ]
        cfg = params.cfg
        cmd += cli_flag(cfg, "min_mapq", "-q")

        with open(str(log[0]), "w") as fh:
            subprocess.check_call(cmd, stdout=fh, stderr=subprocess.STDOUT)


rule cnvkit_reference:
    """Build a pooled reference from normal-sample coverage files."""
    input:
        normal_cnn=lambda wc: [
            config["bb_dir"]
            + f"/{wc.assay_type}/cnvkit/{rep_id}.{rt}coverage.cnn"
            for rep_id, st in zip(
                assay2rep_ids[wc.assay_type],
                assay2sample_types[wc.assay_type],
            ) if st == "normal"
            for rt in ["target", "antitarget"]
        ],
    output:
        reference_cnn=config["bb_dir"] + "/{assay_type}/cnvkit/reference.cnn",
    wildcard_constraints:
        assay_type="bulkWES",
    params:
        cnvkit=_cnvkit_exe,
        reference=config["reference"],
        cfg=_cnvkit_cfg,
    log:
        config["log_dir"] + "/cnvkit_reference.{assay_type}.log",
    conda:
        "../envs/cnvkit.yaml"
    run:
        import subprocess

        cmd = [
            str(params.cnvkit), "reference",
        ]
        cmd += [str(f) for f in input.normal_cnn]
        cmd += ["-f", str(params.reference)]
        cmd += ["-o", str(output.reference_cnn)]
        cfg = params.cfg
        cmd += cli_flag(cfg, "no_gc", "--no-gc", is_bool=True)
        cmd += cli_flag(cfg, "no_edge", "--no-edge", is_bool=True)
        cmd += cli_flag(cfg, "no_rmask", "--no-rmask", is_bool=True)
        cmd += cli_flag(cfg, "male_reference", "--male-reference", is_bool=True)
        cmd += cli_flag(cfg, "cluster", "--cluster", is_bool=True)

        with open(str(log[0]), "w") as fh:
            subprocess.check_call(cmd, stdout=fh, stderr=subprocess.STDOUT)


rule cnvkit_fix:
    """Correct each replicate's coverage against the pooled reference."""
    input:
        target_cnn=config["bb_dir"]
        + "/{assay_type}/cnvkit/{rep_id}.targetcoverage.cnn",
        antitarget_cnn=config["bb_dir"]
        + "/{assay_type}/cnvkit/{rep_id}.antitargetcoverage.cnn",
        reference_cnn=config["bb_dir"] + "/{assay_type}/cnvkit/reference.cnn",
    output:
        cnr=config["bb_dir"] + "/{assay_type}/cnvkit/{rep_id}.cnr",
    wildcard_constraints:
        assay_type="bulkWES",
    params:
        cnvkit=_cnvkit_exe,
        cfg=_cnvkit_cfg,
    log:
        config["log_dir"] + "/cnvkit_fix.{assay_type}_{rep_id}.log",
    conda:
        "../envs/cnvkit.yaml"
    run:
        import subprocess

        cmd = [
            str(params.cnvkit), "fix",
            str(input.target_cnn),
            str(input.antitarget_cnn),
            str(input.reference_cnn),
            "-o", str(output.cnr),
        ]
        cfg = params.cfg
        cmd += cli_flag(cfg, "no_gc", "--no-gc", is_bool=True)
        cmd += cli_flag(cfg, "no_edge", "--no-edge", is_bool=True)
        cmd += cli_flag(cfg, "no_rmask", "--no-rmask", is_bool=True)

        with open(str(log[0]), "w") as fh:
            subprocess.check_call(cmd, stdout=fh, stderr=subprocess.STDOUT)


rule cnvkit_to_window_dp:
    """Convert .cnr files to window.dp.npz + window.tsv.gz for combine_counts."""
    input:
        cnr_files=lambda wc: [
            config["bb_dir"]
            + f"/{wc.assay_type}/cnvkit/{rep_id}.cnr"
            for rep_id in assay2rep_ids[wc.assay_type]
        ],
        sample_file=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/sample_ids.tsv",
        region_bed=config["region_bed"],
        blacklist_bed=config.get("blacklist_bed", []),
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
    log:
        config["log_dir"] + "/cnvkit_to_window_dp.{assay_type}.log",
    conda:
        "../envs/cnvkit.yaml"
    script:
        """../scripts/cnvkit_to_window_dp.py"""
