# Unified Genotyping Pipeline

Unified Genotyping pipeline inputs BAMs from various assay types (bulk WGS/WES, scRNA, scATAC, VISIUM, VISIUM HD 3prime), performs SNP genotyping, phasing, allele counting, bias correction, and segmentation, and outputs segment-level allele count and feature count matrices.

---

## Documentation

| Document | Description |
|----------|-------------|
| [config/config.yaml](config/config.yaml) | Workflow configuration file. |
| [config/samples.tsv](config/samples.tsv) | Template of sample file. |
| [docs/step-by-step.md](docs/step-by-step.md) | Step-by-step guide for running all three workflow modes. |
| [docs/spec.md](docs/spec.md) | Input/output reference. |
| [resources/README.md](resources/README.md) | External data: SNP panels, phasing panels, GTF files, bias BEDs, phasing tools. |
