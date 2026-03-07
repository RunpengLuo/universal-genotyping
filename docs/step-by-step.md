# Step-by-Step Guide

All modes share the same command structure. Fill in paths in `config/config.yaml` before running.

```sh
snakemake --cores <N> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output_dir> \
    --config sample_file=/path/to/samples.tsv sample_id=<PATIENT_ID>
```

`sample_id` selects which patient to process (must match a `SAMPLE` value in the sample sheet). `assay_types` in the config controls which assay types to include.

---

## Mode 1: `bulk_genotyping`

For bulk WGS/WES BAMs. Genotypes SNPs via `bcftools`, phases with Eagle/SHAPEIT/LongPhase, computes allele counts and bias-corrected RDR.

1. Prepare a sample sheet with bulk samples (`assay_type` in `{bulkWGS, bulkWES}`), both normal and tumor.
2. Set `workflow_mode: bulk_genotyping` and `assay_types: ["bulkWGS"]` (or `["bulkWES"]`) in config.
3. **WGS:** Set `window_bed` to a pre-filtered window BED with GC/mappability/replication-timing covariates (e.g., produced by `resources/scripts/build_window_bed.py`). Bias correction uses HMMcopy-style LOWESS via `rd_correct.py`.
4. **WES:** Set `wes_targets_bed` to the capture kit targets BED file and configure `params_cnvkit` in config. The pipeline runs CNVkit (`autobin → coverage → reference → fix`) for read depth normalization, with CNVkit handling GC-content, edge-effect, and repeat-mask corrections internally. The `.cnr` output is then converted to the standard window depth format. Note: `window_bed` and mosdepth are **not** used for WES — set `region_bed` if you need region constraints.
5. Run the pipeline. Primary outputs in `bb_dir/{assay_type}/`: `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele,baf,depth,rdr}.npz`.

---

## Mode 2: `single_cell_genotyping`

For scRNA, scATAC, VISIUM, or VISIUM3prime. Genotypes SNPs from pseudobulk pileup via `cellsnp-lite`, phases, pileups per cell, and bins allele counts.

1. Prepare a sample sheet with non-bulk samples. For multiome, give scRNA and scATAC rows the same `REP_ID`. Fill in `PATH_to_barcodes` and `PATH_to_10x_ranger`.
2. Set `workflow_mode: single_cell_genotyping` and `assay_types` accordingly.
3. Run the pipeline. Primary outputs in `bb_dir/{assay_type}/`: `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele}.npz`, plus `barcodes.tsv.gz` in `allele_dir`.

**Tip:** To skip re-genotyping when matched bulk data was already processed, set `het_snp_vcf` in the config to point to the prior run's `phase/phased_het_snps.vcf.gz`.

---

## Mode 3: `copytyping_preprocess`

For when phased het SNPs and a CNV profile (e.g., from HATCHet3) are already available. Skips genotyping/phasing; aggregates per-cell allele counts onto CNV segments.

1. Prepare a sample sheet with non-bulk samples (same as Mode 2).
2. Set `workflow_mode: copytyping_preprocess`, and provide `het_snp_vcf` and `seg_ucn` in config.
3. Run the pipeline. Primary outputs in `bb_dir/{assay_type}/`: `cnv_segments.tsv`, `{X,Y,D}_count.npz`.
