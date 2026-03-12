# Step-by-Step Guide

## Installation

Requires [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/) >= 7.0.

All tool dependencies are managed via conda environments under `workflow/envs/`:

| Env file | Tools |
|----------|-------|
| `base.yaml` | Python scientific stack (scipy, numpy, pandas, numba, anndata, scanpy, etc.) |
| `tools.yaml` | bcftools, cellsnp-lite, mosdepth, samtools, tabix |
| `phase.yaml` | eagle2, shapeit5, longphase, bcftools, tabix |
| `cnvkit.yaml` | cnvkit (WES only) |
| `snapatac2.yaml` | snapatac2, scanpy (scATAC only) |


```sh
snakemake --profile /path/to/workflow/profile/ \
    --conda-create-envs-only --cores 1 \
    -s /path/to/workflow/Snakefile
```

---

## Running the pipeline

```sh
snakemake --profile /path/to/workflow/profile/ \
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
3. **WGS:** Set `window_bed` to a pre-filtered window BED (e.g., build using `resources/scripts/build_wgs_window_bed.py` or use pre-built `resources/data/windows.1kbp.hg38.bed.gz`).
4. **WES:** Set `wes_targets_bed` and configure `params_cnvkit`. CNVkit handles WES RD preprocessing and GC corrections. `wes_targets_bed` is target capture regions provided by the sequencing platform.
5. Outputs in `bb_dir/{assay_type}/`:
   - `bb.tsv.gz` — bin annotations.
   - `bb.{Tallele,Aallele,Ballele,baf,depth,rdr}.npz` — allele, depth, and RDR matrices.
   - `sample_ids.tsv` — sample metadata.

---

## Mode 2: `single_cell_genotyping`

scRNA, scATAC, VISIUM, or VISIUM3prime. Pseudobulk genotyping via cellsnp-lite, per-cell pileup, binned allele counts.

1. Sample sheet with non-bulk samples. Multiome: same `REP_ID` for scRNA + scATAC. Fill in `PATH_to_barcodes` and `PATH_to_10x_ranger`.
2. Set `workflow_mode: single_cell_genotyping` in config.
3. Outputs in `bb_dir/{assay_type}/`:
   - `bb.tsv.gz` — bin annotations.
   - `bb.{Tallele,Aallele,Ballele,baf}.npz` — allele count and BAF matrices.
   - `multi_snp.tsv.gz` — multi-SNP group annotations.
   - `multi_snp.{Tallele,Aallele,Ballele}.npz` — multi-SNP allele count matrices.
   - `sample_ids.tsv` — sample metadata.

To reuse genotyped and phased SNPs from bulk data, set `het_snp_vcf` to a prior run's `phase/phased_het_snps.vcf.gz`.

---

## Mode 3: `copytyping_preprocess`

Skips genotyping. Aggregates per-cell allele counts onto pre-computed CNV segments.

1. Sample sheet with non-bulk samples.
2. Set `workflow_mode: copytyping_preprocess`. Provide `het_snp_vcf`, `seg_ucn`, `bbc_ucn`, and `bbc_phases` in config.
3. Outputs in `bb_dir/{assay_type}/`:
   - `cnv_segments.tsv` — CNV segment annotations.
   - `{X,Y,D}_count.npz` — per-segment count matrices.
   - `barcodes.tsv.gz` — cell barcode list.
   - `sample_ids.tsv` — sample metadata.

---

## Phasing

Set `phaser` in config to one of:

| Phaser | Config keys needed |
|--------|--------------------|
| `eagle` | `eagle_dir` — path to Eagle2 directory (must contain `tables/` with genetic map files) |
| `shapeit` | `shapeit_dir` — path to SHAPEIT5 directory (must contain `resources/maps/`) |
| `longphase` | `params_longphase` — `min_mapq`, `extra_params` (e.g., `"--pb"` for PacBio) |

Eagle and shapeit require a phasing reference panel (`phasing_panel` in config). Longphase does not use a genetic map or phasing panel.
