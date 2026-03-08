# Installation

## Prerequisites

- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/) (recommended)
- [Snakemake](https://snakemake.readthedocs.io/) >= 7.0

## Setup

Pre-create per-rule conda environments:

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
