##################################################
if workflow_mode == "bulk":
    rule postprocess_bulk:
        input:
            vcfs=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.base.vcf.gz", rep_id=mod2reps["bulkDNA"]),
            tot_mats=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.tag.DP.mtx", rep_id=mod2reps["bulkDNA"]),
            ad_mats=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.tag.AD.mtx", rep_id=mod2reps["bulkDNA"]),
            snp_file=lambda wc: "phase/phased_snps.vcf.gz",
            genome_size=lambda wc: config["genome_size"],
            region_bed=lambda wc: config["region_bed"],
            gmap_file=lambda wc: ("phase/genetic_map.tsv.gz" if require_genetic_map else []),
            block_bed=lambda wc: [] if config.get("block_bed") is None else config["block_bed"],
        output:
            sample_file="allele/bulkDNA/sample_ids.tsv",
            snp_file="allele/bulkDNA/snps.tsv.gz",
            tot_mtx_snp="allele/bulkDNA/snp.Tallele.npz",
            a_mtx_snp="allele/bulkDNA/snp.Aallele.npz",
            b_mtx_snp="allele/bulkDNA/snp.Ballele.npz",
            bb_file="allele/bulkDNA/bb.tsv",
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
            "logs/postprocess.bulkDNA.log"
        script:
            """../scripts/postprocess_bulk.py"""

    rule run_mosdepth_bulk:
        input:
            bam=lambda wc: get_data[("bulkDNA", wc.rep_id)][1],
            bed_file="allele/bulkDNA/bb.bed.gz",
        output:
            reg_file="allele/bulkDNA/out_mosdepth/{rep_id}.regions.bed.gz",
        threads: config["threads"]["mosdepth"]
        params:
            mosdepth=config["mosdepth"],
            sample_name=SAMPLE_ID,
            modality="bulkDNA",
            out_prefix="allele/bulkDNA/out_mosdepth/{rep_id}",
            read_quality=config["params_mosdepth"]["read_quality"],
            extra_params=config["params_mosdepth"].get("extra_params", ""),
        log:
            "logs/mosdepth.bulkDNA_{rep_id}.log"
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
            bb_file="allele/bulkDNA/bb.tsv",
            mos_files=expand("allele/bulkDNA/out_mosdepth/{rep_id}.regions.bed.gz", rep_id=mod2reps.get("bulkDNA", [])),
            reference=config["reference"],
            genome_size=config["genome_size"],
            mappability_file=[] if config.get("mappability_file") is None else config["mappability_file"],
        output:
            dp_mtx_bb="allele/bulkDNA/bb.depth.npz",
            rdr_mtx_bb="allele/bulkDNA/bb.rdr.npz",
        params:
            sample_name=SAMPLE_ID,
            qc_dir="allele/bulkDNA/qc",
            mosdepth_dir="allele/bulkDNA/out_mosdepth",
            gc_correct=config["params_compute_rdr"]["gc_correct"]
        log:
            "logs/compute_rdrs.bulkDNA.log"
        script:
            """../scripts/compute_rdr_bulk.py"""


##################################################
if workflow_mode == "single_cell":
    def single_modality_output(wc):
        has_cn_profile = config.get("seg_ucn") is not None
        data_type = f"{wc.data_type}"
        base_dir = f"allele/{data_type}"
        outputs = {
            "sample_file": base_dir + "/sample_ids.tsv",
            "all_barcodes": base_dir + "/barcodes.txt",
            "qc_dir": directory(base_dir + "/qc"),
            "snp_file": base_dir + "/snps.tsv.gz",
            "tot_mtx_snp": [base_dir + "/snp.Tallele.npz"],
            "a_mtx_snp": [base_dir + "/snp.Aallele.npz"],
            "b_mtx_snp": [base_dir + "/snp.Ballele.npz"],
        }
        if has_cn_profile:
            outputs["cnp_file"] = [base_dir + "/haplotype_blocks.tsv"]
            outputs["y_count"] = [base_dir + "/Y_count.npz"]
            outputs["d_count"] = [base_dir + "/D_count.npz"]
        else:
            outputs["meta_file"] = [base_dir + "/meta.tsv.gz"]
            outputs["tot_mtx_meta"] = [base_dir + "/meta.Tallele.npz"]
            outputs["a_mtx_meta"] = [base_dir + "/meta.Aallele.npz"]
            outputs["b_mtx_meta"] = [base_dir + "/meta.Ballele.npz"]
            outputs["bb_file"] = [base_dir + "/bb.tsv"]
            outputs["tot_mtx_bb"] = [base_dir + "/bb.Tallele.npz"]
            outputs["a_mtx_bb"] = [base_dir + "/bb.Aallele.npz"]
            outputs["b_mtx_bb"] = [base_dir + "/bb.Ballele.npz"]

            # legacy CalicoST inputs
            if data_type in ["VISIUM", "VISIUM3PRIME"]:
                outputs["unique_snp_ids_legacy"] = [base_dir + "/unique_snp_ids.npy"]
                outputs["a_mtx_legacy"] = [base_dir + "/cell_snp_Aallele.npz"]
                outputs["b_mtx_legacy"] = [base_dir + "/cell_snp_Ballele.npz"]
        return outputs

    rule postprocess_single_modality:
        input:
            vcfs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            sample_tsvs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.samples.tsv", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            tot_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            ad_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            block_bed=lambda wc: [] if config.get("block_bed") is None else config["block_bed"],
            seg_ucn = lambda wc: [] if config.get("seg_ucn") is None else config["seg_ucn"],
            bbc_ucn = lambda wc: [] if config.get("bbc_ucn") is None else config["bbc_ucn"],
            bbc_phases = lambda wc: [] if config.get("bbc_phases") is None else config["bbc_phases"],
        output:
            single_modality_output
        wildcard_constraints:
            data_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        params:
            sample_name=SAMPLE_ID,
            modality=lambda wc: wc.data_type,
            data_types=lambda wc: [wc.data_type],
            rep_ids=lambda wc: mod2reps.get(wc.data_type, []),
            min_depth=config["params_postprocess"]["min_depth"],
            gamma=config["params_postprocess"]["gamma"],
            nu=config["params_postprocess"]["nu"],
            min_switchprob=config["params_postprocess"]["min_switchprob"],
            max_switchprob=config["params_postprocess"]["max_switchprob"],
            binom_test=config["params_postprocess"]["binom_test"],
            binom_alpha=config["params_postprocess"]["binom_alpha"],
            binom_margin=config["params_postprocess"]["binom_margin"],
            nsnp_meta=config["params_postprocess"]["nsnp_meta"],
            min_snp_reads=config["params_postprocess"]["min_snp_reads"],
            min_snp_per_block=config["params_postprocess"]["min_snp_per_block"],
        log:
            "logs/postprocess.{data_type}.log"
        script:
            """../scripts/postprocess_single_cell.py"""

    ##################################################
    def multiome_output(wc):
        data_types = ["scRNA", "scATAC"]
        has_cn_profile = config.get("seg_ucn") is not None
        base_dir = f"allele/multiome_{wc.rep_id}"
        outputs = {
            "sample_file": base_dir + "/sample_ids.tsv",
            "all_barcodes": base_dir + "/barcodes.txt",
            "qc_dir": directory(base_dir + "/qc"),
            "snp_file": base_dir + "/snps.tsv.gz",
            "tot_mtx_snp": [base_dir + f"/snp.Tallele.{data_type}.npz" for data_type in data_types],
            "a_mtx_snp": [base_dir + f"/snp.Aallele.{data_type}.npz" for data_type in data_types],
            "b_mtx_snp": [base_dir + f"/snp.Ballele.{data_type}.npz" for data_type in data_types],
        }
        if has_cn_profile:
            outputs["cnp_file"] = [base_dir + "/haplotype_blocks.tsv"]
            outputs["y_count"] = []
            outputs["d_count"] = []
            for data_type in data_types:
                outputs["y_count"].append(base_dir + f"/Y_count.{data_type}.npz")
                outputs["d_count"].append(base_dir + f"/D_count.{data_type}.npz")
        else:
            raise ValueError("adaptive binning for multiome data TODO")
        return outputs
    
    rule postprocess_multiome:
        input:
            vcfs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.base.vcf.gz" for data_type in ["scRNA", "scATAC"]],
            sample_tsvs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.samples.tsv" for data_type in ["scRNA", "scATAC"]],
            tot_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.DP.mtx" for data_type in ["scRNA", "scATAC"]],
            ad_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.AD.mtx" for data_type in ["scRNA", "scATAC"]],
            snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            block_bed=lambda wc: [] if config.get("block_bed") is None else config["block_bed"],
            seg_ucn = lambda wc: [] if config.get("seg_ucn") is None else config["seg_ucn"],
            bbc_ucn = lambda wc: [] if config.get("bbc_ucn") is None else config["bbc_ucn"],
            bbc_phases = lambda wc: [] if config.get("bbc_phases") is None else config["bbc_phases"],
        output:
            multiome_output
        params:
            sample_name=SAMPLE_ID,
            modality="multiome",
            data_types=["scRNA", "scATAC"],
            rep_ids=lambda wc: [wc.rep_id],
            min_depth=config["params_postprocess"]["min_depth"],
            gamma=config["params_postprocess"]["gamma"],
            nu=config["params_postprocess"]["nu"],
            min_switchprob=config["params_postprocess"]["min_switchprob"],
            max_switchprob=config["params_postprocess"]["max_switchprob"],
            binom_test=config["params_postprocess"]["binom_test"],
            binom_alpha=config["params_postprocess"]["binom_alpha"],
            binom_margin=config["params_postprocess"]["binom_margin"],
            nsnp_meta=config["params_postprocess"]["nsnp_meta"],
            min_snp_reads=config["params_postprocess"]["min_snp_reads"],
            min_snp_per_block=config["params_postprocess"]["min_snp_per_block"],
        log:
            "logs/postprocess.multiome.{rep_id}.log"
        script:
            """../scripts/postprocess_single_cell.py"""
