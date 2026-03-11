# Unified Genotyping Pipeline

Unified Genotyping pipeline inputs BAMs from various assay types (bulk WGS/WES, scRNA, scATAC, VISIUM, VISIUM HD 3prime), performs SNP genotyping, phasing, allele counting, bias correction, and segmentation, and outputs segment-level allele count and feature count matrices.

---

## Quick-Start

Requires [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/) and [Snakemake](https://snakemake.readthedocs.io/) >= 7.0.

Modify `profile/config.yaml` for system configuration and conda path, etc.,

```sh
# install necessary conda environment.
snakemake --profile /path/to/profile/ \
    --conda-create-envs-only --cores 1 \
    -s /path/to/workflow/Snakefile

# use `--dry-run` to check input files are formatted properly
snakemake --profile /path/to/profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv sample_id=<PATIENT_ID>
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [config/config.yaml](config/config.yaml) | Workflow configuration file. |
| [config/samples.tsv](config/samples.tsv) | Template of sample file. |
| [docs/step-by-step.md](docs/step-by-step.md) | Step-by-step guide for running all three workflow modes. |
| [docs/spec.md](docs/spec.md) | Input/output reference. |
| [resources/README.md](resources/README.md) | External data: SNP panels, phasing panels, GTF files, bias BEDs, phasing tools. |
