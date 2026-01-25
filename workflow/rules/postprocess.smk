##################################################
rule postprocess_matrix_bulk:
    input:
        vcfs=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.base.vcf.gz", rep_id=mod2reps["bulkDNA"]),
        dp_mats=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.tag.DP.mtx", rep_id=mod2reps["bulkDNA"]),
        ad_mats=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.tag.AD.mtx", rep_id=mod2reps["bulkDNA"]),
        snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
        region_bed=lambda wc: config["region_bed"],
    output:
        sample_file="allele/bulkDNA/samples_ids.tsv",
        info_file="allele/bulkDNA/snp_info.tsv.gz",
        dp_mtx="allele/bulkDNA/snp_matrix.dp.npz",
        alt_mtx="allele/bulkDNA/snp_matrix.alt.npz",
        ref_mtx="allele/bulkDNA/snp_matrix.ref.npz",
    params:
        sample_name=SAMPLE_ID,
        modality="bulkDNA",
        data_types=["bulkDNA"],
        rep_ids=mod2reps.get("bulkDNA", None),
        mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
        min_depth=config["params_postprocess"]["min_depth"],
        gamma=config["params_postprocess"]["gamma"],
    log:
        "logs/postprocess.bulkDNA.log"
    script:
        """../scripts/postprocess_bulk.py"""

##################################################
rule postprocess_matrix_aggregate_one_data_type:
    input:
        vcfs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        sample_tsvs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.samples.tsv", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        dp_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        ad_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
        region_bed=lambda wc: config["region_bed"],
    output:
        sample_file="allele/{data_type}/samples_ids.tsv",
        info_file="allele/{data_type}/snp_info.tsv.gz",
        dp_mtx="allele/{data_type}/snp_matrix.dp.npz",
        alt_mtx="allele/{data_type}/snp_matrix.alt.npz",
        ref_mtx="allele/{data_type}/snp_matrix.ref.npz",
        all_barcodes="allele/{data_type}/barcodes.txt",
        unique_snp_ids="allele/{data_type}/unique_snp_ids.npy",
        a_mtx="allele/{data_type}/cell_snp_Aallele.npz",
        b_mtx="allele/{data_type}/cell_snp_Ballele.npz",
    wildcard_constraints:
        data_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
    params:
        sample_name=SAMPLE_ID,
        modality=lambda wc: wc.data_type,
        data_types=lambda wc: [wc.data_type],
        rep_ids=lambda wc: mod2reps.get(wc.data_type, None),
        mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
    log:
        "logs/postprocess.{data_type}.log"
    script:
        """../scripts/postprocess_nonbulk.py"""


rule postprocess_matrix_multiome_per_replicate:
    input:
        vcfs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.base.vcf.gz" for data_type in ["scRNA", "scATAC"]],
        sample_tsvs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.samples.tsv" for data_type in ["scRNA", "scATAC"]],
        dp_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.DP.mtx" for data_type in ["scRNA", "scATAC"]],
        ad_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.AD.mtx" for data_type in ["scRNA", "scATAC"]],
        snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
        region_bed=lambda wc: config["region_bed"],
    output:
        sample_file="allele/multiome_{rep_id}/samples_ids.tsv",
        all_barcodes="allele/multiome_{rep_id}/barcodes.txt",
        info_file="allele/multiome_{rep_id}/snp_info.tsv.gz",
        unique_snp_ids="allele/multiome_{rep_id}/unique_snp_ids.npy",
        dp_mtx=["allele/multiome_{rep_id}/" + f"snp_matrix.{data_type}.dp.npz" for data_type in ["scRNA", "scATAC"]],
        alt_mtx=["allele/multiome_{rep_id}/" + f"snp_matrix.{data_type}.alt.npz" for data_type in ["scRNA", "scATAC"]],
        ref_mtx=["allele/multiome_{rep_id}/" + f"snp_matrix.{data_type}.ref.npz" for data_type in ["scRNA", "scATAC"]],
        a_mtx=["allele/multiome_{rep_id}/" + f"cell_snp_Aallele.{data_type}.npz" for data_type in ["scRNA", "scATAC"]],
        b_mtx=["allele/multiome_{rep_id}/" + f"cell_snp_Ballele.{data_type}.npz" for data_type in ["scRNA", "scATAC"]],
    params:
        sample_name=SAMPLE_ID,
        modality="multiome",
        data_types=["scRNA", "scATAC"],
        rep_ids=lambda wc: [wc.rep_id],
        mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
    log:
        "logs/postprocess.multiome_{rep_id}.log"
    script:
        """../scripts/postprocess_nonbulk.py"""
