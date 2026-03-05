##################################################
# RDR estimation rules (bulk only)
# Controlled by config["rdr_method"]: "bin" or "window"
##################################################

if workflow_mode == "bulk_genotyping":

    _rdr_method = config.get("rdr_method", "bin")
    assert _rdr_method in ("bin", "window"), (
        f"rdr_method must be 'bin' or 'window', got '{_rdr_method}'"
    )

    ##################################################
    # Bin-level RDR (adaptive bins from SNP pipeline)
    # Active when rdr_method == "bin"
    ##################################################

    if _rdr_method == "bin":

        _rdr_cfg_bin = config["params_compute_rdr"]
        assert _rdr_cfg_bin.get("bias_bed") is not None, (
            "rdr_method='bin' requires params_compute_rdr.bias_bed to be set"
        )

        rule run_mosdepth_bulk:
            input:
                bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
                bed_file=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/raw/bb.raw.bed.gz",
            output:
                mosdepth_file=config["bb_dir"]
                + "/{assay_type}/out_mosdepth/{rep_id}.regions.bed.gz",
            threads: config["threads"]["mosdepth"]
            wildcard_constraints:
                assay_type="(bulkDNA|bulkWGS|bulkWES)",
            params:
                sample_name=SAMPLE_ID,
                assay_type=lambda wc: wc.assay_type,
                mosdepth=config["mosdepth"],
                out_prefix=config["bb_dir"] + "/{assay_type}/out_mosdepth/{rep_id}",
                read_quality=config["params_mosdepth"]["read_quality"],
                extra_params=config["params_mosdepth"].get("extra_params", ""),
            log:
                config["log_dir"] + "/run_mosdepth_bulk/run_mosdepth_bulk.{assay_type}_{rep_id}.log",
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
                sample_file=lambda wc: config["bb_dir"]
                + f"/{wc.assay_type}/sample_ids.tsv",
                bb_file=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/raw/bb.raw.tsv.gz",
                mosdepth_files=lambda wc: branch(
                    assay2rep_ids.get(wc.assay_type) is not None,
                    then=[
                        config["bb_dir"]
                        + f"/{wc.assay_type}/out_mosdepth/{rep_id}.regions.bed.gz"
                        for rep_id in assay2rep_ids[wc.assay_type]
                    ],
                    otherwise=[],
                ),
                baf_mtx_bb=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/raw/bb.raw.baf.npz",
                bias_bed=_rdr_cfg_bin["bias_bed"],
                genome_size=config["genome_size"],
                region_bed=branch(
                    config.get("region_bed") is None,
                    then=[],
                    otherwise=config["region_bed"],
                ),
            output:
                dp_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.depth.npz",
                rdr_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.rdr.npz",
                corr_factors=config["bb_dir"] + "/{assay_type}/raw/bb.raw.corr_factors.tsv.gz",
                qc_dir=directory(config["qc_dir"] + "/{assay_type}/compute_rdr_bulk/"),
            wildcard_constraints:
                assay_type="(bulkDNA|bulkWGS|bulkWES)",
            params:
                sample_name=SAMPLE_ID,
                mosdepth_dir=lambda wc: config["bb_dir"]
                + f"/{wc.assay_type}/out_mosdepth",
                correction_method=_rdr_cfg_bin.get("correction_method", "quantile"),
                gc_correct=_rdr_cfg_bin.get("gc_correct", True),
                rt_correct=_rdr_cfg_bin.get("rt_correct", False),
                chromosomes=config["chromosomes"],
            log:
                config["log_dir"] + "/compute_rdr_bulk.{assay_type}.log",
            script:
                """../scripts/compute_rdr_bulk_bin.py"""

    ##################################################
    # Window-based RDR (fixed-width windows + HMMcopy loess)
    # Active when rdr_method == "window"
    ##################################################

    if _rdr_method == "window":

        _rdr_cfg = config["params_compute_rdr"]
        assert _rdr_cfg.get("bias_bed") is not None, (
            "rdr_method='window' requires params_compute_rdr.bias_bed to be set"
        )

        rule run_mosdepth_bulk_window:
            input:
                bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
            output:
                mosdepth_file=config["bb_dir"]
                + "/{assay_type}/out_mosdepth_window/{rep_id}.regions.bed.gz",
            threads: config["threads"]["mosdepth"]
            wildcard_constraints:
                assay_type="(bulkDNA|bulkWGS|bulkWES)",
            params:
                mosdepth=config["mosdepth"],
                out_prefix=config["bb_dir"] + "/{assay_type}/out_mosdepth_window/{rep_id}",
                read_quality=config["params_mosdepth"]["read_quality"],
                extra_params=config["params_mosdepth"].get("extra_params", ""),
                window_size=_rdr_cfg.get("window_size", 1000),  # default 1000 bp
            log:
                config["log_dir"]
                + "/run_mosdepth_bulk_window/run_mosdepth_bulk_window.{assay_type}_{rep_id}.log",
            shell:
                r"""
                {params.mosdepth} \
                    -t {threads} \
                    -Q {params.read_quality} \
                    --by {params.window_size} \
                    {params.extra_params} \
                    "{params.out_prefix}" \
                    "{input.bam}" > {log} 2>&1
                cat "{params.out_prefix}.mosdepth.summary.txt" >> {log} 2>/dev/null || true
                """

        rule compute_rdr_bulk_window:
            input:
                sample_file=lambda wc: config["bb_dir"]
                + f"/{wc.assay_type}/sample_ids.tsv",
                bb_file=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/raw/bb.raw.tsv.gz",
                bias_bed=_rdr_cfg["bias_bed"],
                mosdepth_files=lambda wc: [
                    config["bb_dir"]
                    + f"/{wc.assay_type}/out_mosdepth_window/{rep_id}.regions.bed.gz"
                    for rep_id in assay2rep_ids[wc.assay_type]
                ],
                genome_size=config["genome_size"],
                region_bed=lambda wc: config["region_bed"],
                blacklist_bed=branch(
                    config.get("blacklist_bed") is None,
                    then=[],
                    otherwise=config["blacklist_bed"],
                ),
            output:
                dp_mtx=config["bb_dir"] + "/{assay_type}/window.depth.npz",
                rdr_mtx=config["bb_dir"] + "/{assay_type}/window.rdr.npz",
                rdr_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.rdr.npz",
                dp_mtx_bb=config["bb_dir"] + "/{assay_type}/raw/bb.raw.depth.npz",
                bins_tsv=config["bb_dir"] + "/{assay_type}/window.bins.tsv.gz",
                corr_factors=config["bb_dir"] + "/{assay_type}/raw/bb.raw.corr_factors.tsv.gz",
                qc_dir=directory(config["qc_dir"] + "/{assay_type}/compute_rdr_bulk/"),
            wildcard_constraints:
                assay_type="(bulkDNA|bulkWGS|bulkWES)",
            params:
                sample_name=SAMPLE_ID,
                mosdepth_dir=lambda wc: config["bb_dir"]
                + f"/{wc.assay_type}/out_mosdepth_window",
                chromosomes=config["chromosomes"],
                samplesize=_rdr_cfg.get("samplesize", 50000),
                routlier=_rdr_cfg.get("routlier", 0.01),
                doutlier=_rdr_cfg.get("doutlier", 0.001),
                min_mappability=_rdr_cfg.get("min_mappability", 0.9),
                gc_correct=_rdr_cfg.get("gc_correct", True),
                rt_correct=_rdr_cfg.get("rt_correct", False),
            log:
                config["log_dir"] + "/compute_rdr_bulk_window.{assay_type}.log",
            script:
                """../scripts/compute_rdr_bulk_window.py"""
