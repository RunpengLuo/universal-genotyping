# Input / Output Reference

## Sample File

TSV listing all data from one patient. Template at `config/samples.tsv`.

| Column | Required | Description |
|--------|----------|-------------|
| `SAMPLE` | Yes | Patient ID. Set `sample_id` in config to select which patient to process. |
| `REP_ID` | Yes | Replicate ID. Unique per row, except multiome pairs (scRNA + scATAC share a `REP_ID`). |
| `assay_type` | Yes | One of: `bulkWGS`, `bulkWES`, `scATAC`, `scRNA`, `VISIUM`, `VISIUM3prime`. |
| `sample_type` | Yes | `normal` or `tumor`. |
| `PATH_to_bam` | Yes | Path to `.bam` file. |
| `PATH_to_barcodes` | Non-bulk | Path to 10x `barcodes.tsv.gz`. |
| `PATH_to_10x_ranger` | Non-bulk | Path to 10x Cell Ranger / Space Ranger `outs/` directory. |
| `PATH_to_ref_annotations` | Optional | Reference annotation TSV with columns `BARCODE` and `<ref_label>`. |

---

## Output Files

Output directories are set in `config.yaml`: `snp_dir`, `phase_dir`, `pileup_dir`, `allele_dir`, `bb_dir`, `qc_dir`, `log_dir`.

### 1. SNP Genotyping (`snp_dir/`)

| File | Description |
|------|-------------|
| `chr{chrname}.vcf.gz` | Per-chromosome VCF of bi-allelic het/hom-alt SNPs via `bcftools` (bulk) or `cellsnp-lite` pseudobulk (single-cell). |
| `pseudobulk_{modality}/cellSNP.*` | Pseudobulk pileup outputs from `cellsnp-lite`. *(single-cell only)* |

### 2. SNP Phasing (`phase_dir/`)

| File | Description |
|------|-------------|
| `chr{chrname}.vcf.gz` | Per-chromosome phased VCF from the configured phaser. |
| `phased_het_snps.vcf.gz(.tbi)` | Genome-wide phased het SNPs (concatenated); input for downstream pileup. |
| `genetic_map.tsv.gz` | Parsed genetic map (cM positions per SNP) for switch-prob estimation. *(eagle/shapeit only)* |

### 3. Per-Replicate Pileup (`pileup_dir/`)

One subdirectory per `{assay_type}_{rep_id}` containing `cellsnp-lite` output: `cellSNP.base.vcf.gz`, `cellSNP.tag.{DP,AD}.mtx`, `cellSNP.samples.tsv`.

### 4. Allele Counts (`allele_dir/{assay_type}/`)

Intermediate allele-count matrices before binning.

| File | Description |
|------|-------------|
| `snps.tsv.gz` | SNP annotation table (chr, pos, ref/alt, phase, block, feature). |
| `snp.{Tallele,Aallele,Ballele}.npz` | Sparse matrices (SNPs x samples/cells): total, A-allele, B-allele counts. |
| `sample_ids.tsv` | Sample metadata (rep IDs, sample types). |
| `barcodes.tsv.gz` | Cell barcode / sample ID list matching matrix columns. |
| `{assay_type}.h5ad` | AnnData with cells x bins count matrix from gene expression or ATAC tiles. *(single-cell only)* |

### 5. Final Bins (`bb_dir/{assay_type}/`)

Primary pipeline outputs.

**Bulk (`bulk_genotyping`):**

| File | Description |
|------|-------------|
| `bb.tsv.gz` | Genomic bin annotation (chr, start, end, SNP count, features). |
| `bb.{Tallele,Aallele,Ballele,baf}.npz` | Allele count and BAF matrices (bins x samples). |
| `bb.depth.npz` | Read depth matrix (bins x samples), aggregated from corrected windows. |
| `bb.rdr.npz` | Read depth ratio matrix (tumor/normal or median-centered). |
| `bb.bed.gz` | BED3 file of bin intervals. |
| `sample_ids.tsv` | Sample metadata. |

Window-level intermediates:

| File | Description |
|------|-------------|
| `window.dp.npz` | Bias-corrected per-window depth matrix. |
| `window.df.tsv.gz` | Window annotation with GC/mappability/replication-timing columns. |

**Single-cell (`single_cell_genotyping`):**

| File | Description |
|------|-------------|
| `bb.tsv.gz` | Genomic bin annotation. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele count matrices (bins x cells). |

**CNV segments (`copytyping_preprocess`):**

| File | Description |
|------|-------------|
| `cnv_segments.tsv` | CNV segment table with copy numbers from input segmentation. |
| `{X,Y,D}_count.npz` | Per-cell x per-segment: total, A-allele, and depth count matrices. |

### 6. QC (`qc_dir/`)

QC plots are produced during `rd_correct` (bias correction QC), `combine_counts` (bulk), and `combine_counts_nonbulk` (single-cell). Output directories: `qc_dir/{assay_type}/rd_correction/` (bulk bias correction), `qc_dir/{assay_type}/segmentation/` (bulk), and `qc_dir/{assay_type}/combine_counts/` (single-cell).

---

## Key TSV Columns

### `snps.tsv.gz`

`#CHR`, `POS`, `POS0`, `START`, `END`, `GT`, `PHASE` (0 = B-allele is ALT, 1 = B-allele is REF), `region_id` (string, `CHR:START-END` of enclosing region), `feature_id`, `feature_type` (exon/intron/intergenic).

### `bb.tsv.gz`

`#CHR`, `START`, `END`, `START0`, `END0`, `switchprobs`, `region_id`, `PS`, `#SNPS`, `BLOCKSIZE`, `feature_id`, `feature_type`, `bb_id`.

### `sample_ids.tsv`

`SAMPLE`, `SAMPLE_NAME`, `REP_ID`, `sample_type`.

### `barcodes.tsv.gz`

Single column, no header. Each row: `{BARCODE}_{REP_ID}`.
