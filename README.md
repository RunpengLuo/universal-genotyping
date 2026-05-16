# Universal Genotyping Pipeline (pre-release)

Universal Genotyping Pipeline (pre-release) inputs BAMs from various assay types (bulk LR/SR WGS/WES, scRNA-seq, scATAC-seq, VISIUM, VISIUM HD 3prime), performs panel-based SNP genotyping, population-based/read-based haplotype phasing, read counting, sequencing bias correction, and segmentation, and outputs genomic segment-level count matrices.

---

## Documentation

| Document | Description |
|----------|-------------|
| [config/config.yaml](config/config.yaml) | Default configuration (auto-loaded by Snakefile). |
| [resources/templates/config.yaml](resources/templates/config.yaml) | Minimal user config template (run-specific overrides). |
| [resources/templates/samples.tsv](resources/templates/samples.tsv) | Template sample sheet. |
| [docs/README.md](docs/README.md) | Step-by-step guide for running all three workflow modes. |
| [docs/tutorials/bulk_genotyping.hg38_wgs.md](docs/tutorials/bulk_genotyping.hg38_wgs.md) | Tutorial: bulk WGS genotyping on hg38. |
| [docs/reference.md](docs/reference.md) | Input/output reference. |
| [resources/README.md](resources/README.md) | External resources. |
