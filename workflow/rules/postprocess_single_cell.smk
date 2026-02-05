##################################################
if workflow_mode == "single_cell":

    rule postprocess_single_modality:
        input:
            vcfs=lambda wc: expand(
                "pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz",
                data_type=wc.data_type,
                rep_id=mod2reps[wc.data_type],
            ),
            sample_tsvs=lambda wc: expand(
                "pileup/{data_type}_{rep_id}/cellSNP.samples.tsv",
                data_type=wc.data_type,
                rep_id=mod2reps[wc.data_type],
            ),
            tot_mats=lambda wc: expand(
                "pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx",
                data_type=wc.data_type,
                rep_id=mod2reps[wc.data_type],
            ),
            ad_mats=lambda wc: expand(
                "pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx",
                data_type=wc.data_type,
                rep_id=mod2reps[wc.data_type],
            ),
            snp_file=lambda wc: branch(
                run_genotype_snps,
                then="phase/phased_snps.vcf.gz",
                otherwise=config["ref_snp_file"],
            ),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            gmap_file=lambda wc: (
                "phase/genetic_map.tsv.gz" if require_genetic_map else []
            ),
            block_bed=lambda wc: branch(
                config.get("block_bed") is None, then=[], otherwise=config["block_bed"]
            ),
            seg_ucn=lambda wc: branch(
                config.get("seg_ucn") is None, then=[], otherwise=config["seg_ucn"]
            ),
            bbc_ucn=lambda wc: branch(
                config.get("bbc_ucn") is None, then=[], otherwise=config["bbc_ucn"]
            ),
            bbc_phases=lambda wc: branch(
                config.get("bbc_phases") is None,
                then=[],
                otherwise=config["bbc_phases"],
            ),
        output:
            sample_file="allele/{data_type}/sample_ids.tsv",
            all_barcodes="allele/{data_type}/barcodes.txt",
            qc_dir=directory("allele/{data_type}/qc"),
            snp_file="allele/{data_type}/snps.tsv.gz",
            tot_mtx_snp=["allele/{data_type}/snp.Tallele.npz"],
            a_mtx_snp=["allele/{data_type}/snp.Aallele.npz"],
            b_mtx_snp=["allele/{data_type}/snp.Ballele.npz"],
        wildcard_constraints:
            data_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        params:
            sample_name=SAMPLE_ID,
            modality=lambda wc: wc.data_type,
            data_types=lambda wc: [wc.data_type],
            rep_ids=lambda wc: mod2reps.get(wc.data_type, []),
            nu=config["params_postprocess"]["nu"],
            min_switchprob=config["params_postprocess"]["min_switchprob"],
            max_switchprob=config["params_postprocess"]["max_switchprob"],
            binom_test=config["params_postprocess"]["binom_test"],
            binom_alpha=config["params_postprocess"]["binom_alpha"],
            binom_margin=config["params_postprocess"]["binom_margin"],
            nsnp_meta=config["params_postprocess"]["nsnp_meta"],
            min_snp_reads=config["params_postprocess"]["min_snp_reads"],
            min_snp_per_block=config["params_postprocess"]["min_snp_per_block"],
            # optional outputs, when CNP is provided
            has_cn_profile=lambda wc: config.get("seg_ucn") is not None,
            cnp_file=lambda wc: f"allele/{wc.data_type}/haplotype_blocks.tsv",
            y_count=lambda wc: [f"allele/{wc.data_type}/Y_count.npz"],
            d_count=lambda wc: [f"allele/{wc.data_type}/D_count.npz"],
            # otherwise, meta-SNP outputs
            meta_file=lambda wc: f"allele/{wc.data_type}/meta.tsv.gz",
            tot_mtx_meta=lambda wc: [f"allele/{wc.data_type}/meta.Tallele.npz"],
            a_mtx_meta=lambda wc: [f"allele/{wc.data_type}/meta.Aallele.npz"],
            b_mtx_meta=lambda wc: [f"allele/{wc.data_type}/meta.Ballele.npz"],
            # and BB segment outputs
            bb_file=lambda wc: f"allele/{wc.data_type}/bb.tsv.gz",
            tot_mtx_bb=lambda wc: [f"allele/{wc.data_type}/bb.Tallele.npz"],
            a_mtx_bb=lambda wc: [f"allele/{wc.data_type}/bb.Aallele.npz"],
            b_mtx_bb=lambda wc: [f"allele/{wc.data_type}/bb.Ballele.npz"],
            # legacy CalicoST snp-level data
            unique_snp_ids_legacy=lambda wc: f"allele/{wc.data_type}/unique_snp_ids.npy",
            a_mtx_legacy=lambda wc: f"allele/{wc.data_type}/cell_snp_Aallele.npz",
            b_mtx_legacy=lambda wc: f"allele/{wc.data_type}/cell_snp_Ballele.npz",
        log:
            "logs/postprocess.{data_type}.log",
        script:
            """../scripts/postprocess_single_cell.py"""

    ##################################################
    rule postprocess_multiome:
        input:
            vcfs=lambda wc: [
                f"pileup/{data_type}_{wc.rep_id}/cellSNP.base.vcf.gz"
                for data_type in ["scRNA", "scATAC"]
            ],
            sample_tsvs=lambda wc: [
                f"pileup/{data_type}_{wc.rep_id}/cellSNP.samples.tsv"
                for data_type in ["scRNA", "scATAC"]
            ],
            tot_mats=lambda wc: [
                f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.DP.mtx"
                for data_type in ["scRNA", "scATAC"]
            ],
            ad_mats=lambda wc: [
                f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.AD.mtx"
                for data_type in ["scRNA", "scATAC"]
            ],
            snp_file=lambda wc: branch(
                run_genotype_snps,
                then="phase/phased_snps.vcf.gz",
                otherwise=config["ref_snp_file"],
            ),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            gmap_file=lambda wc: (
                "phase/genetic_map.tsv.gz" if require_genetic_map else []
            ),
            block_bed=lambda wc: branch(
                config.get("block_bed") is None, then=[], otherwise=config["block_bed"]
            ),
            seg_ucn=lambda wc: branch(
                config.get("seg_ucn") is None, then=[], otherwise=config["seg_ucn"]
            ),
            bbc_ucn=lambda wc: branch(
                config.get("bbc_ucn") is None, then=[], otherwise=config["bbc_ucn"]
            ),
            bbc_phases=lambda wc: branch(
                config.get("bbc_phases") is None,
                then=[],
                otherwise=config["bbc_phases"],
            ),
        output:
            sample_file="allele/multiome_{rep_id}/sample_ids.tsv",
            all_barcodes="allele/multiome_{rep_id}/barcodes.txt",
            qc_dir=directory("allele/multiome_{rep_id}/qc"),
            snp_file="allele/multiome_{rep_id}/snps.tsv.gz",
            tot_mtx_snp=expand(
                "allele/multiome_{{rep_id}}/snp.Tallele.{data_type}.npz",
                data_type=["scRNA", "scATAC"],
            ),
            a_mtx_snp=expand(
                "allele/multiome_{{rep_id}}/snp.Aallele.{data_type}.npz",
                data_type=["scRNA", "scATAC"],
            ),
            b_mtx_snp=expand(
                "allele/multiome_{{rep_id}}/snp.Ballele.{data_type}.npz",
                data_type=["scRNA", "scATAC"],
            ),
        params:
            sample_name=SAMPLE_ID,
            modality="multiome",
            data_types=["scRNA", "scATAC"],
            rep_ids=lambda wc: [wc.rep_id],
            nu=config["params_postprocess"]["nu"],
            min_switchprob=config["params_postprocess"]["min_switchprob"],
            max_switchprob=config["params_postprocess"]["max_switchprob"],
            binom_test=config["params_postprocess"]["binom_test"],
            binom_alpha=config["params_postprocess"]["binom_alpha"],
            binom_margin=config["params_postprocess"]["binom_margin"],
            nsnp_meta=config["params_postprocess"]["nsnp_meta"],
            min_snp_reads=config["params_postprocess"]["min_snp_reads"],
            min_snp_per_block=config["params_postprocess"]["min_snp_per_block"],
            # optional outputs, when CNP is provided
            has_cn_profile=lambda wc: config.get("seg_ucn") is not None,
            cnp_file=lambda wc: "allele/multiome_{wc.rep_id}/haplotype_blocks.tsv",
            y_count=lambda wc: [
                f"allele/multiome_{wc.rep_id}/Y_count.{data_type}.npz"
                for data_type in ["scRNA", "scATAC"]
            ],
            d_count=lambda wc: [
                f"allele/multiome_{wc.rep_id}/D_count.{data_type}.npz"
                for data_type in ["scRNA", "scATAC"]
            ],
            # otherwise, meta-SNP outputs TODO
            # meta_file=lambda wc: "allele/multiome_{wc.rep_id}/meta.tsv.gz",
            # tot_mtx_meta=lambda wc: "allele/multiome_{wc.rep_id}/meta.Tallele.npz",
            # a_mtx_meta=lambda wc: "allele/multiome_{wc.rep_id}/meta.Aallele.npz",
            # b_mtx_meta=lambda wc: "allele/multiome_{wc.rep_id}/meta.Ballele.npz",
            # # and BB segment outputs
            # bb_file=lambda wc: f"allele/multiome_{wc.rep_id}/bb.tsv.gz",
            # tot_mtx_bb=lambda wc: f"allele/multiome_{wc.rep_id}/bb.Tallele.npz",
            # a_mtx_bb=lambda wc: "allele/multiome_{wc.rep_id}/bb.Aallele.npz",
            # b_mtx_bb=lambda wc: "allele/multiome_{wc.rep_id}/bb.Ballele.npz",
        log:
            "logs/postprocess.multiome.{rep_id}.log",
        script:
            """../scripts/postprocess_single_cell.py"""
