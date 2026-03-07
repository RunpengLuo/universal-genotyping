# Installation

## Prerequisites

- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/) (recommended)
- [Snakemake](https://snakemake.readthedocs.io/) >= 7.0

## Per-rule conda environments (recommended)

The pipeline ships per-rule environment files under `workflow/envs/`, avoiding dependency conflicts between packages (e.g., cnvkit and snapatac2).

| Env file | Key packages | Used by |
|---|---|---|
| `workflow/envs/base.yaml` | python, scipy, numpy, pandas, numba, scanpy, squidpy, statsmodels | Most Python script rules |
| `workflow/envs/tools.yaml` | bcftools, samtools, cellsnp-lite, mosdepth | Shell-based bioinformatics rules |
| `workflow/envs/cnvkit.yaml` | cnvkit, numpy, pandas | CNVkit WES pipeline rules |
| `workflow/envs/snapatac2.yaml` | snapatac2, scanpy, anndata | `process_atac_fragments` |

### Snakemake profile

A local profile at `profile/config.yaml` sets default CLI flags so you don't have to repeat them:

```yaml
use-conda: true          # enable per-rule conda envs
conda-prefix: .snakemake/conda   # shared env cache (relative to repo root)
conda-frontend: mamba    # faster env creation
```

Just add `--profile profile/` to your snakemake command. To override `conda-prefix` (e.g., for a shared filesystem), pass `--conda-prefix /absolute/path` — CLI flags always take precedence.

### Pre-create environments

To create all environments upfront without running the pipeline:

```sh
snakemake --profile profile/ --conda-create-envs-only --cores 1 \
    -s workflow/Snakefile
```

No `--configfile` is needed — the Snakefile detects the absence of config and exposes dummy rules referencing all four env files so Snakemake can discover and create them.

By default, environments are installed under the path set by `conda-prefix` in `profile/config.yaml` (`.snakemake/conda/`, relative to the working directory). To use a custom location (e.g., a shared filesystem on a cluster), pass `--conda-prefix` on the command line:

```sh
snakemake --profile profile/ --conda-create-envs-only --cores 1 \
    -s workflow/Snakefile \
    --conda-prefix /path/to/shared/conda-envs
```

CLI flags always take precedence over the profile.

## Manual single-environment install (development)

For development, install everything into one environment via `environment.yaml`:

```sh
conda env create -f environment.yaml
conda activate genotyping-env
```

Then run without the profile. Note: may have dependency conflicts; per-rule approach is preferred for production.
