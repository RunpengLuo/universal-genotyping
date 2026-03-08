# Unified Genotyping Pipeline

A Snakemake pipeline for unified SNP genotyping, phasing, and allele counting across multiple assay types (bulk WGS/WES, scRNA, scATAC, VISIUM, VISIUM HD 3prime) from the same patient. Supports bulk genotyping, single-cell genotyping, and copy-typing preprocessing modes.

## Table of Contents

- [Installation](#installation)
- [Quick-Start](#quick-start)
- [Workflow Overview](#workflow-overview)
- [Documentation](#documentation)

---

## Installation

Requires [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/) and [Snakemake](https://snakemake.readthedocs.io/) >= 7.0.

Setup environments, modify `profile/config.yaml` if required, all 
```sh
snakemake --profile profile/ --conda-create-envs-only --cores 1 \
    -s workflow/Snakefile --conda-prefix /absolute/path
```

| Env file | Packages | Rules |
|---|---|---|
| `workflow/envs/base.yaml` | python, scipy, numpy, pandas, numba, scanpy, statsmodels | Python scripts |
| `workflow/envs/tools.yaml` | bcftools, samtools, cellsnp-lite, mosdepth | Bioinformatics tools |
| `workflow/envs/cnvkit.yaml` | cnvkit | WES depth correction |
| `workflow/envs/snapatac2.yaml` | snapatac2, scanpy | ATAC fragment processing |

---

## Quick-Start

```sh
# dry-run to check input files are formatted properly
snakemake --cores 1 \
    --dry-run \
    --profile profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv sample_id=<PATIENT_ID>

# run the pipeline
snakemake --cores <num_cores> \
    --profile profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv sample_id=<PATIENT_ID>
```

---

## Workflow Overview

### Bulk Genotyping

![Bulk genotyping rulegraph](docs/rulegraph.bulk_genotyping.png)

### Single-Cell Genotyping

![Single-cell genotyping rulegraph](docs/rulegraph.single_cell_genotyping.png)

---

## Documentation

| Document | Description |
|----------|-------------|
| [config/config.yaml](config/config.yaml) | Workflow configuration file. |
| [config/samples.tsv](config/samples.tsv) | Template of sample file. |
| [docs/step-by-step.md](docs/step-by-step.md) | Step-by-step guide for running all three workflow modes. |
| [docs/spec.md](docs/spec.md) | Input/output reference. |
| [resources/README.md](resources/README.md) | External data: SNP panels, phasing panels, GTF files, bias BEDs, phasing tools. |
