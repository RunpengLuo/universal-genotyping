# Tutorial: Bulk WGS Genotyping (hg38)

This tutorial walks through running the `bulk_genotyping` pipeline on paired normal/tumor bulk WGS samples aligned to GRCh38.

## 1. Prerequisites

Install [conda](https://github.com/conda-forge/miniforge) and [Snakemake](https://snakemake.readthedocs.io/) (>= 7.0). We recommend using the [libmamba solver](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) for faster dependency resolution (`conda config --set solver libmamba`).

### External data

| Resource | Description |
|----------|-------------|
| **SNP panel + targets + phasing panel** | 1000 Genomes phase 3 VCF — we recommend the [1kGP n=3,202 high-coverage](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) panel |
| **Eagle2** | Phasing tool ([download](https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz)) — includes genetic maps in `tables/` |
| **SHAPEIT5** | Phasing tool ([GitHub](https://github.com/odelaneau/shapeit5)) — includes genetic maps in `resources/maps/` |
| **GTF annotation** | GENCODE v38 ([download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)) |

Download and prepare the 1kGP n=3,202 panel (requires bcftools, bgzip, tabix):

```bash
bash resources/scripts/process_1kGP_3202_panel.sh /path/to/1kGP_3202
```

This downloads the raw VCFs and produces three directories under the output path:

| Output | Config key | Contents |
|--------|------------|----------|
| `snps.vcf.gz` | `snp_panel` | Merged SNP-only VCF |
| `target_positions/` | `snp_targets` | `target.chr{1..22,X}.pos.gz` + `.tbi` |
| `phasing_panel/` | `phasing_panel` | `chr{1..22,X}.genotypes.bcf` + `.csi` |

### Pre-built files

The following are already included in `resources/data/` and do not need to be downloaded:

- `windows.1kbp.hg38.bed.gz` — 1 kb window BED with GC, mappability, and replication timing
- `hg38-blacklist.v2.bed.gz` — ENCODE blacklist v2
- `hg38.regions.bed` — chromosome region definitions

## 2. Sample sheet

Create a TSV file (e.g., `samples.tsv`) with one row per BAM. For bulk WGS you need five columns:

```tsv
SAMPLE	REP_ID	assay_type	sample_type	PATH_to_bam
HT001	N1	bulkWGS	normal	/data/HT001/normal.bam
HT001	T1	bulkWGS	tumor	/data/HT001/tumor.bam
```

- `SAMPLE` — patient identifier (must match `sample_id` in config)
- `REP_ID` — unique per row
- `assay_type` — `bulkWGS`
- `sample_type` — `normal` or `tumor`
- `PATH_to_bam` — absolute path to the BAM file (must be indexed)

## 3. Config (`config/config.yaml`)

Copy the template config and edit the sections below. Only the fields that need changing from defaults are shown.

### Workflow mode

```yaml
workflow_mode: "bulk_genotyping"
assay_types: ["bulkWGS"]
```

### Input

```yaml
sample_id: HT001
sample_file: /path/to/samples.tsv
chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
```

### Reference

```yaml
reference_version: hg38
reference: /path/to/reference.fasta
genome_size: resources/data/hg38.regions.bed
region_bed: resources/data/hg38.regions.bed
window_bed: resources/data/windows.1kbp.hg38.bed.gz
blacklist_bed: resources/data/hg38-blacklist.v2.bed.gz
gtf_file: /path/to/gencode.v38.annotation.gtf.gz
```

### Population panels

```yaml
snp_panel: /path/to/snp_panel.vcf.gz
snp_targets: /path/to/snp_targets
phasing_panel: /path/to/phasing_panel
```

The `phasing_panel` directory should contain per-chromosome BCF files.

### Phasing

Pick one of Eagle2 or SHAPEIT5. Both require a genetic map (gmap) bundled with the tool distribution.

**Eagle2** — set `eagle_dir` to the extracted Eagle2 directory. The pipeline looks for the gmap at `{eagle_dir}/tables/genetic_map_hg38_withX.txt.gz`.

```yaml
phaser: "eagle"
eagle_dir: /path/to/Eagle_v2.4.1
```

**SHAPEIT5** — set `shapeit_dir` to the SHAPEIT5 repository root. The pipeline looks for per-chromosome gmaps at `{shapeit_dir}/resources/maps/b38/chr{chrname}.b38.gmap.gz`.

```yaml
phaser: "shapeit"
shapeit_dir: /path/to/shapeit5
```

### Count reads
For bulk WGS data, by default we apply quadratic median regression to correct GC content and replication timing biases while being robust to copy-number signals. Refer to `<qc_dir>/<assay_type>/rd_correction/*`, check read-depth plots and KDE scatter plot for diagonstic purpose.

### Combine counts

For bulk WGS, `min_snp_reads` and `min_snp_per_block` balanced the bin size and phase-switch errors from reference phasing, we recommend to leave them as default, but if the number of bins are pretty small due to low coverage, you may decrease them accordingly.

```yaml
params_combine_counts:
  min_snp_reads: 1000
  min_snp_per_block: 10
```

## 4. Profile and running

### Create the Snakemake profile

The repository includes a profile at `profile/config.yaml`. You may ddjust `cores` to match your machine.

```yaml
cores: 16
use-conda: true
conda-prefix: .snakemake/conda
printshellcmds: true
show-failed-logs: true
verbose: true
rerun-triggers:
  - mtime
  - input
```

### Pre-create conda environments

```bash
snakemake --profile profile/ \
    --conda-create-envs-only --cores 1 \
    -s workflow/Snakefile \
    --configfile config/config.yaml
```

### Run the pipeline

```bash
snakemake --profile profile/ \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --directory /path/to/output_dir \
    --config sample_file=/path/to/samples.tsv sample_id=HT001
```

The `--config` flags on the command line override values in `config.yaml`, so you can reuse the same config for different patients.

## 5. Results

Final outputs required by HATCHet3 can be found in `<out_dir>/<bb_dir>/bulkWGS/`. See [docs/spec.md](../spec.md) for the full output specification including intermediate files.
