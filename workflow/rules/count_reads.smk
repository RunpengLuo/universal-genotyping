##################################################
# Mosdepth on pre-filtered windows + per-window RD bias correction (bulk only)
# window_bed (from build_window_bed.py) contains pre-filtered windows with
# GC/MAP/REPLI/region_id/is_target columns. A headerless 3-column BED
# is derived for mosdepth --by.
##################################################

if workflow_mode == "bulk_genotyping":

    _rdr_cfg = config["params_count_reads"]

    # -----------------------------------------------------------------
    # Derive headerless windows BED from window_bed for mosdepth
    # -----------------------------------------------------------------

    rule window_bed_to_3bed:
        """Extract headerless 3-column BED from window_bed for mosdepth --by."""
        input:
            window_bed=config["window_bed"],
        output:
            mosdepth_bed=config["bb_dir"] + "/windows.bed.gz",
        log:
            config["log_dir"] + "/window_bed_to_3bed.log",
        run:
            import gzip, logging
            import pandas as pd
            logging.basicConfig(filename=str(log[0]), level=logging.INFO,
                                format="%(asctime)s %(levelname)s %(message)s")

            df = pd.read_table(str(input.window_bed), sep="\t", usecols=["#CHR", "START", "END"])
            logging.info(f"Extracted {len(df)} windows from {input.window_bed}")
            with gzip.open(str(output.mosdepth_bed), "wt") as fh:
                df.to_csv(fh, sep="\t", header=False, index=False)

    # -----------------------------------------------------------------
    # Mosdepth
    # -----------------------------------------------------------------

    rule run_mosdepth:
        input:
            bam=lambda wc: get_data[(wc.assay_type, wc.rep_id)][1],
            windows_bed=config["bb_dir"] + "/windows.bed.gz",
        output:
            mosdepth_file=config["bb_dir"]
            + "/{assay_type}/out_mosdepth/{rep_id}.regions.bed.gz",
        threads: config["threads"]["mosdepth"]
        wildcard_constraints:
            assay_type="(bulkWGS|bulkWES)",
        params:
            mosdepth=config["mosdepth"],
            out_prefix=config["bb_dir"] + "/{assay_type}/out_mosdepth/{rep_id}",
            read_quality=config["params_mosdepth"]["read_quality"],
            extra_params=config["params_mosdepth"].get("extra_params", ""),
        log:
            config["log_dir"]
            + "/run_mosdepth/run_mosdepth.{assay_type}_{rep_id}.log",
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

    # -----------------------------------------------------------------
    # Per-window RD bias correction (HMMcopy-style LOWESS)
    # -----------------------------------------------------------------

    rule rd_correct:
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
        output:
            dp_corrected=config["bb_dir"] + "/{assay_type}/window.depth_corrected.npz",
            window_df=config["bb_dir"] + "/{assay_type}/window.df.tsv.gz",
            qc_dir=directory(config["qc_dir"] + "/{assay_type}/rd_correction/"),
        wildcard_constraints:
            assay_type="(bulkWGS|bulkWES)",
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
        script:
            """../scripts/rd_correct.py"""
