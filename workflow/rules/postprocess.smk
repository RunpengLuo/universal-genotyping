##################################################
if workflow_mode == "bulk":
    rule postprocess_matrix_bulk:
        input:
            vcfs=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.base.vcf.gz", rep_id=mod2reps["bulkDNA"]),
            tot_mats=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.tag.DP.mtx", rep_id=mod2reps["bulkDNA"]),
            ad_mats=lambda wc: expand("pileup/bulkDNA_{rep_id}/cellSNP.tag.AD.mtx", rep_id=mod2reps["bulkDNA"]),
            snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
        output:
            sample_file="allele/bulkDNA/samples_ids.tsv",
            info_file="allele/bulkDNA/snp_info.tsv.gz",
            bed_file="allele/bulkDNA/snp_positions.bed.gz",
            tot_mtx="allele/bulkDNA/snp_matrix.tot.npz",
            alt_mtx="allele/bulkDNA/snp_matrix.alt.npz",
            ref_mtx="allele/bulkDNA/snp_matrix.ref.npz",
            a_mtx="allele/bulkDNA/snp_matrix.A.npz",
            b_mtx="allele/bulkDNA/snp_matrix.B.npz",
        params:
            sample_name=SAMPLE_ID,
            modality="bulkDNA",
            data_types=["bulkDNA"],
            rep_ids=mod2reps.get("bulkDNA", []),
            mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
            min_depth=config["params_postprocess"]["min_depth"],
            gamma=config["params_postprocess"]["gamma"],
        log:
            "logs/postprocess.bulkDNA.log"
        script:
            """../scripts/postprocess_bulk.py"""

    rule run_mosdepth_bulk:
        input:
            bam=lambda wc: get_data[("bulkDNA", wc.rep_id)][1],
            bed_file="allele/bulkDNA/snp_positions.bed.gz",
        output:
            reg_file="allele/bulkDNA/out_mosdepth/{rep_id}.regions.bed.gz",
        threads: config["threads"]["mosdepth"]
        params:
            mosdepth=config["mosdepth"],
            sample_name=SAMPLE_ID,
            modality="bulkDNA",
            out_prefix="allele/bulkDNA/out_mosdepth/{rep_id}",
            read_quality=config["params_mosdepth"].get("read_quality", 11),
            extra_params=config["params_mosdepth"].get("extra_params", ""),
        log:
            "logs/mosdepth_counting.bulkDNA_{rep_id}.log"
        shell:
            r"""
            {params.mosdepth} \
                -t {threads} \
                -Q {params.read_quality} \
                --by "{input.bed_file}" \
                {params.extra_params} \
                "{params.out_prefix}" \
                "{input.bam}"
            """

    rule count_reads_bulk:
        input:
            sample_file="allele/bulkDNA/samples_ids.tsv",
            info_file="allele/bulkDNA/snp_info.tsv.gz",
            mos_files=expand("allele/bulkDNA/out_mosdepth/{rep_id}.regions.bed.gz", rep_id=mod2reps.get("bulkDNA", [])),
        output:
            dp_mtx="allele/bulkDNA/snp_matrix.dp.npz",
        params:
            mosdepth_dir="allele/bulkDNA/out_mosdepth",
        log:
            "logs/count_reads.bulkDNA.log"
        run:
            import pandas as pd
            import numpy as np

            print("merge mosdepth results across multi-samples")
            samples_df = pd.read_table(input.sample_file, sep="\t")
            rep_ids = samples_df["REP_ID"].astype(str).tolist()
            snp_info = pd.read_table(input.info_file, sep="\t")
            for rep_id in rep_ids:
                mos_file = f"{params.mosdepth_dir}/{rep_id}.regions.bed.gz"
                mos_df = pd.read_table(
                    mos_file, sep="\t", header=None, names=["#CHR", "START", "END", rep_id]
                )
                dup = mos_df.duplicated(subset=["#CHR", "START", "END"], keep=False)
                assert not dup.any(), f"{mos_file} has duplicate intervals: {int(dup.sum())}"
                assert len(mos_df) == len(snp_info), f"{mos_file} rowcount={len(mos_df)} != snp_info rowcount={len(snp_info)}"
                snp_info = snp_info.merge(
                    mos_df, how="left", on=["#CHR", "START", "END"], sort=False
                )
            dp_mat = snp_info.loc[rep_ids].to_numpy(dtype=np.float32)
            np.savez_compressed(output.dp_mtx, mat=dp_mat)


##################################################
if workflow_mode == "single_cell":
    rule postprocess_matrix_aggregate_one_data_type:
        input:
            vcfs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.base.vcf.gz", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            sample_tsvs=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.samples.tsv", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            tot_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.DP.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            ad_mats=lambda wc: expand("pileup/{data_type}_{rep_id}/cellSNP.tag.AD.mtx", data_type=wc.data_type, rep_id=mod2reps[wc.data_type]),
            snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
        output:
            sample_file="allele/{data_type}/samples_ids.tsv",
            info_file="allele/{data_type}/snp_info.tsv.gz",
            bed_file="allele/{data_type}/snp_positions.bed.gz",
            all_barcodes="allele/{data_type}/barcodes.txt",
            tot_mtx="allele/{data_type}/snp_matrix.tot.npz",
            alt_mtx="allele/{data_type}/snp_matrix.alt.npz",
            ref_mtx="allele/{data_type}/snp_matrix.ref.npz",
            a_mtx="allele/{data_type}/snp_matrix.A.npz",
            b_mtx="allele/{data_type}/snp_matrix.B.npz",
            # legacy CalicoST outputs
            unique_snp_ids_legacy="allele/{data_type}/unique_snp_ids.npy",
            a_mtx_legacy="allele/{data_type}/cell_snp_Aallele.npz",
            b_mtx_legacy="allele/{data_type}/cell_snp_Ballele.npz",
        wildcard_constraints:
            data_type="(scRNA|scATAC|VISIUM|VISIUM3prime)",
        params:
            sample_name=SAMPLE_ID,
            modality=lambda wc: wc.data_type,
            data_types=lambda wc: [wc.data_type],
            rep_ids=lambda wc: mod2reps.get(wc.data_type, []),
            mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
        log:
            "logs/postprocess.{data_type}.log"
        script:
            """../scripts/postprocess_nonbulk.py"""

    ##################################################
    rule postprocess_matrix_multiome_per_replicate:
        input:
            vcfs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.base.vcf.gz" for data_type in ["scRNA", "scATAC"]],
            sample_tsvs=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.samples.tsv" for data_type in ["scRNA", "scATAC"]],
            tot_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.DP.mtx" for data_type in ["scRNA", "scATAC"]],
            ad_mats=lambda wc: [f"pileup/{data_type}_{wc.rep_id}/cellSNP.tag.AD.mtx" for data_type in ["scRNA", "scATAC"]],
            snp_file=lambda wc: ("phase/phased_snps.vcf.gz" if run_genotype_snps else config["ref_snp_file"]),
            region_bed=lambda wc: config["region_bed"],
            genome_size=lambda wc: config["genome_size"],
        output:
            sample_file="allele/multiome_{rep_id}/samples_ids.tsv",
            info_file="allele/multiome_{rep_id}/snp_info.tsv.gz",
            bed_file="allele/multiome_{rep_id}/snp_positions.bed.gz",
            all_barcodes="allele/multiome_{rep_id}/barcodes.txt",
            tot_mtx=["allele/multiome_{rep_id}/" + f"{data_type}_snp_matrix.tot.npz" for data_type in ["scRNA", "scATAC"]],
            alt_mtx=["allele/multiome_{rep_id}/" + f"{data_type}_snp_matrix.alt.npz" for data_type in ["scRNA", "scATAC"]],
            ref_mtx=["allele/multiome_{rep_id}/" + f"{data_type}_snp_matrix.ref.npz" for data_type in ["scRNA", "scATAC"]],
            a_mtx=["allele/multiome_{rep_id}/" + f"{data_type}_snp_matrix.A.npz" for data_type in ["scRNA", "scATAC"]],
            b_mtx=["allele/multiome_{rep_id}/" + f"{data_type}_snp_matrix.B.npz" for data_type in ["scRNA", "scATAC"]],
        params:
            sample_name=SAMPLE_ID,
            modality="multiome",
            data_types=["scRNA", "scATAC"],
            rep_ids=lambda wc: [wc.rep_id],
            mask_out_of_region=config["params_postprocess"]["mask_out_of_region"],
        log:
            "logs/postprocess.multiome.{rep_id}.log"
        script:
            """../scripts/postprocess_nonbulk.py"""
