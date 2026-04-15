# Unified Genotyping Pipeline (Under development)

Unified Genotyping pipeline inputs BAMs from various assay types (bulk WGS/WES, scRNA, scATAC, VISIUM, VISIUM HD 3prime), performs SNP genotyping, phasing, allele counting, bias correction, and segmentation, and outputs segment-level allele count and feature count matrices.

---

## Documentation

| Document | Description |
|----------|-------------|
| [config/config.yaml](config/config.yaml) | Default configuration (auto-loaded by Snakefile). |
| [templates/config.yaml](templates/config.yaml) | Minimal user config template (run-specific overrides). |
| [templates/samples.tsv](templates/samples.tsv) | Template sample sheet. |
| [docs/README.md](docs/README.md) | Step-by-step guide for running all three workflow modes. |
| [docs/tutorials/bulk_genotyping.hg38_wgs.md](docs/tutorials/bulk_genotyping.hg38_wgs.md) | Tutorial: bulk WGS genotyping on hg38. |
| [docs/reference.md](docs/reference.md) | Input/output reference. |
| [resources/README.md](resources/README.md) | External resources. |
