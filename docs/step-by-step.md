# Step-by-Step Guide

## Running the pipeline

```sh
snakemake --cores <N> \
    --profile profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output_dir> \
    --config sample_file=/path/to/samples.tsv sample_id=<PATIENT_ID>
```

`sample_id` must match a `SAMPLE` value in the sample sheet. `--profile profile/` enables per-rule conda envs.

---

## Mode 1: `bulk_genotyping`

Bulk WGS/WES. Genotypes SNPs (bcftools), phases, computes allele counts and bias-corrected RDR.

1. Sample sheet with `assay_type` = `bulkWGS` or `bulkWES`, including normal and tumor.
2. Set `workflow_mode: bulk_genotyping` in config.
3. **WGS:** Set `window_bed` to a pre-filtered window BED (e.g., from `resources/scripts/build_window_bed.py`).
4. **WES:** Set `wes_targets_bed` and configure `params_cnvkit`. CNVkit handles all bias correction internally; `window_bed`/mosdepth are not used.
5. Outputs in `bb_dir/{assay_type}/`: `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele,baf,depth,rdr}.npz`.

---

## Mode 2: `single_cell_genotyping`

scRNA, scATAC, VISIUM, or VISIUM3prime. Pseudobulk genotyping via cellsnp-lite, per-cell pileup, binned allele counts.

1. Sample sheet with non-bulk samples. Multiome: same `REP_ID` for scRNA + scATAC. Fill in `PATH_to_barcodes` and `PATH_to_10x_ranger`.
2. Set `workflow_mode: single_cell_genotyping` in config.
3. Outputs in `bb_dir/{assay_type}/`: `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele}.npz`.

To reuse bulk genotyping results, set `het_snp_vcf` to a prior run's `phase/phased_het_snps.vcf.gz`.

---

## Mode 3: `copytyping_preprocess`

Skips genotyping. Aggregates per-cell allele counts onto pre-computed CNV segments.

1. Sample sheet with non-bulk samples.
2. Set `workflow_mode: copytyping_preprocess`. Provide `het_snp_vcf` and `seg_ucn` in config.
3. Outputs in `bb_dir/{assay_type}/`: `cnv_segments.tsv`, `{X,Y,D}_count.npz`.
