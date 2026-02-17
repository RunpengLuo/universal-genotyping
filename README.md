# Unified Genotyping Pipeline

A Snakemake pipeline for unified SNP genotyping, phasing, and allele counting across multiple assay types (bulk WGS/WES, scRNA, scATAC, VISIUM, VISIUM HD 3prime) from the same patient. Supports bulk genotyping, single-cell genotyping, and copy-typing preprocessing modes.

## Table of Contents

- [Quick-Start](#quick-start)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Sample File Format](#sample-file-format)
- [Documentation](#documentation)

---

## Quick-Start

```sh
mamba env create -f ./environment.yaml -p /path/to/envs/genotyping_env
conda activate /path/to/envs/genotyping_env

# dry-run to check input files are formatted properly
snakemake --cores 1 \
    --dry-run \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv

# run the pipeline
snakemake --cores <num_cores> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv

# generate report
snakemake --cores <num_cores> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv \
    --report report.html
```

---

## Dependencies

Refer to `environment.yaml` for the complete list of Python/conda dependencies. Additional phasing tools must be installed separately:

1. **Eagle2** — https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
2. **SHAPEIT5** — https://github.com/odelaneau/shapeit5
3. **LongPhase** — https://github.com/twolinin/longphase

---

## Usage

See [docs/step-by-step.md](docs/step-by-step.md) for a step-by-step guide to running each of the three workflow modes. To configure the pipeline, edit `config/config.yaml` following the inline comments.

---

## Sample File Format

The sample file is a TSV listing all data from the same patient. A template is provided at `config/samples.tsv`.

| Column | Required | Explanation |
|--------|----------|-------------|
| `SAMPLE` | Yes | Patient ID. All rows must have the same value. |
| `REP_ID` | Yes | Replicate ID. Must be unique per row, except for multiome pairs (scRNA + scATAC share a `REP_ID`). |
| `assay_type` | Yes | One of: `bulkWGS`, `bulkWES`, `scATAC`, `scRNA`, `VISIUM`, `VISIUM3prime`. |
| `sample_type` | Yes | `normal` or `tumor`. |
| `PATH_to_bam` | Yes | Path to `.bam` file. |
| `PATH_to_barcodes` | Non-bulk | Path to 10x `barcodes.tsv.gz`. Found under `outs/filtered_feature_bc_matrix/`. |
| `PATH_to_10x_ranger` | Non-bulk | Path to 10x Cell Ranger / Space Ranger `outs/` directory. Must contain `filtered_feature_bc_matrix*` and `{atac_}*fragments.tsv.gz` (scATAC). |
| `PATH_to_ref_annotations` | Optional | Reference annotation TSV with columns `BARCODE` and `<ref_label>`. |

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/step-by-step.md](docs/step-by-step.md) | Step-by-step guide for running all three workflow modes. |
| [docs/output.md](docs/output.md) | Complete reference for all output files and their columns. |
| [docs/bulk_mode.md](docs/bulk_mode.md) | Detailed notes on the bulk genotyping workflow. |
| [docs/single_cell_mode.md](docs/single_cell_mode.md) | Detailed notes on the single-cell genotyping workflow. |
| [docs/resources.md](docs/resources.md) | External data resources (SNP panels, phasing panels, GTF files). |
| [docs/terminology.md](docs/terminology.md) | Definitions for key terms and cellsnp-lite output fields. |
