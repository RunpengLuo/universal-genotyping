##################################################
rule postprocess_matrix_bulk:
    input:
        vcfs=lambda wc: expand("pileup/bulk_DNA_{rep_id}/cellSNP.base.vcf.gz", rep_id=mod2reps["bulk_DNA"]),
        dp_mats=lambda wc: expand("pileup/bulk_DNA_{rep_id}/cellSNP.tag.DP.mtx", rep_id=mod2reps["bulk_DNA"]),
        ad_mats=lambda wc: expand("pileup/bulk_DNA_{rep_id}/cellSNP.tag.AD.mtx", rep_id=mod2reps["bulk_DNA"]),
        snp_file=lambda wc: ("phase/phased_het_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
        region_bed=lambda wc: config["region_bed"],
    output:
        sample_file="allele/bulk_DNA/samples_ids.tsv",
        info_file="allele/bulk_DNA/snp_info.tsv.gz",
        dp_mtx="allele/bulk_DNA/snp_matrix.dp.npz",
        alt_mtx="allele/bulk_DNA/snp_matrix.alt.npz",
        ref_mtx="allele/bulk_DNA/snp_matrix.ref.npz",
    params:
        sample_name=SAMPLE_ID,
        modality="bulk_DNA",
        data_types=["bulk_DNA"],
        rep_ids=mod2reps.get("bulk_DNA", None),
        mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
        min_depth=config["params_postprocess"]["min_depth"],
        gamma=config["params_postprocess"]["gamma"],
    script:
        """../scripts/postprocess_bulk.py"""

##################################################
rule postprocess_matrix_aggregate_one_data_type:
    input:
        vcfs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        sample_tsvs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.samples.tsv", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        dp_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        ad_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
        snp_file=lambda wc: ("phase/phased_het_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
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
        data_type="^(sc_GEX|sc_ATAC|VISIUM|VISIUM_3PRIME)$",
    params:
        sample_name=SAMPLE_ID,
        modality=lambda wc: wc.data_type,
        data_types=lambda wc: [wc.data_type],
        rep_ids=lambda wc: mod2reps.get(wc.data_type, None),
        mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
    script:
        """../scripts/postprocess_nonbulk.py"""


rule postprocess_matrix_multiome_per_replicate:
    input:
        vcfs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.base.vcf.gz" for data_type in ["sc_GEX", "sc_ATAC"]],
        sample_tsvs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.samples.tsv" for data_type in ["sc_GEX", "sc_ATAC"]],
        dp_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.DP.mtx" for data_type in ["sc_GEX", "sc_ATAC"]],
        ad_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.AD.mtx" for data_type in ["sc_GEX", "sc_ATAC"]],
        snp_file=lambda wc: ("phase/phased_het_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
        region_bed=lambda wc: config["region_bed"],
    output:
        sample_file="allele/multiome_{rep_id}/samples_ids.tsv",
        all_barcodes="allele/multiome_{rep_id}/barcodes.txt",
        info_file="allele/multiome_{rep_id}/snp_info.tsv.gz",
        unique_snp_ids="allele/multiome_{rep_id}/unique_snp_ids.npy",
        dp_mtx=["allele/multiome_{rep_id}/" + f"snp_matrix.{data_type}.dp.npz" for data_type in ["sc_GEX", "sc_ATAC"]],
        alt_mtx=["allele/multiome_{rep_id}/" + f"snp_matrix.{data_type}.alt.npz" for data_type in ["sc_GEX", "sc_ATAC"]],
        ref_mtx=["allele/multiome_{rep_id}/" + f"snp_matrix.{data_type}.ref.npz" for data_type in ["sc_GEX", "sc_ATAC"]],
        a_mtx=["allele/multiome_{rep_id}/" + f"cell_snp_Aallele.{data_type}.npz" for data_type in ["sc_GEX", "sc_ATAC"]],
        b_mtx=["allele/multiome_{rep_id}/" + f"cell_snp_Ballele.{data_type}.npz" for data_type in ["sc_GEX", "sc_ATAC"]],
    params:
        sample_name=SAMPLE_ID,
        modality="multiome",
        data_types=["sc_GEX", "sc_ATAC"],
        rep_ids=lambda wc: [wc.rep_id],
        mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
    script:
        """../scripts/postprocess_nonbulk.py"""
