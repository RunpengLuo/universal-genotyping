##################################################
if workflow_mode == "bulk_genotyping":

    rule phase_and_concat_bulk:
        input:
            vcfs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.base.vcf.gz"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            sample_tsvs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.samples.tsv"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            tot_mtxs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.tag.DP.mtx"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            ad_mtxs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.tag.AD.mtx"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            snp_vcf=lambda wc: branch(
                run_genotype_snps,
                then=config["phase_dir"] + "/phased_het_snps.vcf.gz",
                otherwise=config["het_snp_vcf"],
            ),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            gtf_file=lambda wc: config["gtf_file"],
            blacklist_bed=lambda wc: branch(
                config.get("blacklist_bed") is None,
                then=[],
                otherwise=config["blacklist_bed"],
            ),
        output:
            all_barcodes=config["allele_dir"] + "/{assay_type}/barcodes.tsv.gz",
            snp_info=config["allele_dir"] + "/{assay_type}/snps.tsv.gz",
            tot_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Tallele.npz",
            a_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Aallele.npz",
            b_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Ballele.npz",
            sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
            qc_dir=directory(config["qc_dir"] + "/{assay_type}/phase_and_concat/"),
        wildcard_constraints:
            assay_type="(bulkDNA|bulkWGS|bulkWES)",
        params:
            sample_name=SAMPLE_ID,
            assay_type=lambda wc: wc.assay_type,
            rep_ids=lambda wc: assay2rep_ids[wc.assay_type],
            sample_types=lambda wc: assay2sample_types[wc.assay_type],
            min_depth=config["params_postprocess"]["min_depth"],
            gamma=config["params_postprocess"]["gamma"],
            exon_only=config["params_postprocess"]["exon_only"],
        log:
            config["log_dir"] + "/phase_and_concat.{assay_type}.log",
        script:
            """../scripts/phase_and_concat.py"""


##################################################
if workflow_mode in ["single_cell_genotyping", "copytyping_preprocess"]:

    rule process_rna_anndata:
        input:
            barcodes=lambda wc: [
                get_data[(wc.assay_type, rid)][0]
                for rid in assay2rep_ids[wc.assay_type]
            ],
            ranger_dirs=lambda wc: [
                get_data[(wc.assay_type, rid)][2]
                for rid in assay2rep_ids[wc.assay_type]
            ],
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            gene_blacklist_file=lambda wc: branch(
                config.get("gene_blacklist_file") is None,
                then=[],
                otherwise=config["gene_blacklist_file"],
            ),
            gtf_file=lambda wc: config["gtf_file"],
        output:
            h5ad_file=config["bb_dir"] + "/{assay_type}/{assay_type}.h5ad",
            cell_types=config["bb_dir"] + "/{assay_type}/cell_types.tsv.gz",
        params:
            assay_type=lambda wc: wc.assay_type,
            rep_ids=lambda wc: assay2rep_ids[wc.assay_type],
            sample_types=lambda wc: assay2sample_types[wc.assay_type],
            min_frac_barcodes=lambda wc: config["params_postprocess"][
                "min_frac_barcodes"
            ],
            gene_id_colname=lambda wc: config["params_postprocess"]["gene_id_colname"],
            rep2ref_annotation=lambda wc: rep2ref_annotation,
            ref_label=lambda wc: config["ref_label"],
        wildcard_constraints:
            assay_type="(scRNA|VISIUM|VISIUM3prime)",
        log:
            config["log_dir"] + "/process_rna_anndata.{assay_type}.log",
        script:
            """../scripts/process_rna_anndata.py"""

    rule process_atac_fragments:
        input:
            barcodes=lambda wc: [
                get_data[(wc.assay_type, rid)][0]
                for rid in assay2rep_ids[wc.assay_type]
            ],
            ranger_dirs=lambda wc: [
                get_data[(wc.assay_type, rid)][2]
                for rid in assay2rep_ids[wc.assay_type]
            ],
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            gene_blacklist_file=lambda wc: branch(
                config.get("gene_blacklist_file") is None,
                then=[],
                otherwise=config["gene_blacklist_file"],
            ),
            gtf_file=lambda wc: config["gtf_file"],
        output:
            h5ad_file=config["bb_dir"] + "/{assay_type}/{assay_type}.h5ad",
            cell_types=config["bb_dir"] + "/{assay_type}/cell_types.tsv.gz",
        params:
            assay_type=lambda wc: wc.assay_type,
            rep_ids=lambda wc: assay2rep_ids[wc.assay_type],
            sample_types=lambda wc: assay2sample_types[wc.assay_type],
            tile_width=lambda wc: config["params_postprocess"]["tile_width"],
            rep2ref_annotation=lambda wc: rep2ref_annotation,
            ref_label=lambda wc: config["ref_label"],
        wildcard_constraints:
            assay_type="scATAC",
        log:
            config["log_dir"] + "/process_atac_fragments.{assay_type}.log",
        script:
            """../scripts/process_atac_fragments.py"""

    rule phase_and_concat_single_cell:
        input:
            vcfs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.base.vcf.gz"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            sample_tsvs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.samples.tsv"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            tot_mtxs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.tag.DP.mtx"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            ad_mtxs=lambda wc: [
                config["pileup_dir"] + f"/{wc.assay_type}_{rep_id}/cellSNP.tag.AD.mtx"
                for rep_id in assay2rep_ids[wc.assay_type]
            ],
            snp_vcf=lambda wc: branch(
                run_genotype_snps,
                then=config["phase_dir"] + "/phased_het_snps.vcf.gz",
                otherwise=config["het_snp_vcf"],
            ),
            h5ad_file=lambda wc: config["bb_dir"]
            + "/{assay_type}/{assay_type}.h5ad",
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
            gtf_file=lambda wc: config["gtf_file"],
            blacklist_bed=lambda wc: branch(
                config.get("blacklist_bed") is None,
                then=[],
                otherwise=config["blacklist_bed"],
            ),
        output:
            all_barcodes=config["allele_dir"] + "/{assay_type}/barcodes.tsv.gz",
            snp_info=config["allele_dir"] + "/{assay_type}/snps.tsv.gz",
            tot_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Tallele.npz",
            a_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Aallele.npz",
            b_mtx_snp=config["allele_dir"] + "/{assay_type}/snp.Ballele.npz",
            sample_file=config["bb_dir"] + "/{assay_type}/sample_ids.tsv",
            qc_dir=directory(config["qc_dir"] + "/{assay_type}/phase_and_concat/"),
            unique_snp_ids=config["allele_dir"] + "/{assay_type}/unique_snp_ids.npy",
            cell_snp_Aallele=config["allele_dir"] + "/{assay_type}/cell_snp_Aallele.npz",
            cell_snp_Ballele=config["allele_dir"] + "/{assay_type}/cell_snp_Ballele.npz",
        wildcard_constraints:
            assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        params:
            sample_name=SAMPLE_ID,
            assay_type=lambda wc: wc.assay_type,
            rep_ids=lambda wc: assay2rep_ids[wc.assay_type],
            sample_types=lambda wc: assay2sample_types[wc.assay_type],
            min_depth=config["params_postprocess"]["min_depth"],
            gamma=config["params_postprocess"]["gamma"],
            exon_only=config["params_postprocess"]["exon_only"],
        log:
            config["log_dir"] + "/phase_and_concat.{assay_type}.log",
        script:
            """../scripts/phase_and_concat.py"""


rule adaptive_binning:
    input:
        snp_info=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snps.tsv.gz",
        tot_mtx_snp=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/snp.Tallele.npz",
        a_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Aallele.npz",
        b_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Ballele.npz",
        sample_file=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
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
        bed_file=config["bb_dir"] + "/{assay_type}/bb.bed.gz",
        qc_dir=directory(config["qc_dir"] + "/{assay_type}/adaptive_binning/"),
    params:
        sample_name=SAMPLE_ID,
        assay_type=lambda wc: wc.assay_type,
        nu=config["params_postprocess"]["nu"],
        min_switchprob=config["params_postprocess"]["min_switchprob"],
        max_switchprob=config["params_postprocess"]["max_switchprob"],
        switchprob_ps=config["params_postprocess"]["switchprob_ps"],
        binom_test=config["params_postprocess"]["binom_test"],
        binom_alpha=config["params_postprocess"]["binom_alpha"],
        binom_margin=config["params_postprocess"]["binom_margin"],
        nsnp_multi=config["params_postprocess"]["nsnp_multi"],
        min_snp_reads=config["params_postprocess"]["min_snp_reads"],
        min_snp_per_block=config["params_postprocess"]["min_snp_per_block"],
    threads: 1
    log:
        config["log_dir"] + "/adaptive_binning.{assay_type}.log",
    script:
        """../scripts/adaptive_binning.py"""


rule cnv_segmentation:
    input:
        snp_info=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snps.tsv.gz",
        tot_mtx_snp=lambda wc: config["allele_dir"]
        + f"/{wc.assay_type}/snp.Tallele.npz",
        a_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Aallele.npz",
        b_mtx_snp=lambda wc: config["allele_dir"] + f"/{wc.assay_type}/snp.Ballele.npz",
        sample_file=lambda wc: config["bb_dir"] + f"/{wc.assay_type}/sample_ids.tsv",
        all_barcodes=config["allele_dir"] + "/{assay_type}/barcodes.tsv.gz",
        h5ad_file=config["bb_dir"] + "/{assay_type}/{assay_type}.h5ad",
        region_bed=lambda wc: config["region_bed"],
        genome_size=lambda wc: config["genome_size"],
        gtf_file=lambda wc: config["gtf_file"],
        seg_ucn=lambda wc: config["seg_ucn"],
        bbc_ucn=lambda wc: config["bbc_ucn"],
        bbc_phases=lambda wc: config["bbc_phases"],
    output:
        cnv_segments=config["bb_dir"] + "/{assay_type}/cnv_segments.tsv",
        x_count=config["bb_dir"] + "/{assay_type}/X_count.npz",
        y_count=config["bb_dir"] + "/{assay_type}/Y_count.npz",
        d_count=config["bb_dir"] + "/{assay_type}/D_count.npz",
        barcodes_out=config["bb_dir"] + "/{assay_type}/barcodes.tsv.gz",
        qc_dir=directory(config["qc_dir"] + "/{assay_type}/cnv_segmentation/"),
    wildcard_constraints:
        assay_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    params:
        sample_name=SAMPLE_ID,
        assay_type=lambda wc: wc.assay_type,
        feature_type=lambda wc: assay_type2feature_type[wc.assay_type],
    threads: 1
    log:
        config["log_dir"] + "/cnv_segmentation.{assay_type}.log",
    script:
        """../scripts/cnv_segmentation.py"""
