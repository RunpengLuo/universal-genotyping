# Output Files

Files are organized by workflow stage. Paths use the directory prefixes defined in `config.yaml` (`snp_dir`, `phase_dir`, `pileup_dir`, `allele_dir`, `qc_dir`). `{assay_type}` is one of `bulkDNA`, `bulkWGS`, `bulkWES`, `scATAC`, `scRNA`, `VISIUM`, or `VISIUM3prime`; `{rep_id}` is the replicate identifier from the sample sheet.

---

## 1. SNP Genotyping (`snp_dir/`)

| File | Description |
|------|-------------|
| `chr{chrname}.vcf.gz` | Per-chromosome VCF of bi-allelic het and hom-alt SNPs called from the BAM(s) via `bcftools` (bulk) or `cellsnp-lite` pseudobulk pileup (single-cell). |
| `pseudobulk_{modality}/cellSNP.base.vcf.gz` | Pseudobulk SNP VCF produced by `cellsnp-lite` across all BAMs for a given modality (DNA or RNA). *(single-cell modes only)* |
| `pseudobulk_{modality}/cellSNP.samples.tsv` | Sample list associated with the pseudobulk genotyping run. |
| `pseudobulk_{modality}/cellSNP.tag.DP.mtx` | Total depth sparse matrix from the pseudobulk pileup. |
| `pseudobulk_{modality}/cellSNP.tag.AD.mtx` | Alt-allele depth sparse matrix from the pseudobulk pileup. |

---

## 2. SNP Phasing (`phase_dir/`)

| File | Description |
|------|-------------|
| `chr{chrname}.vcf.gz` | Per-chromosome phased VCF produced by the selected phaser (SHAPEIT4, Eagle, or LongPhase). |
| `phased_het_snps.vcf.gz` | Genome-wide VCF of phased heterozygous SNPs (genotypes `0\|1` or `1\|0`) concatenated across all chromosomes; used as the SNP reference for downstream pileup. |
| `phased_het_snps.vcf.gz.tbi` | Tabix index for `phased_het_snps.vcf.gz`. |
| `genetic_map.tsv.gz` | Parsed genetic map table (cM positions per SNP) used to compute haplotype-switch probabilities in adaptive binning. *(eagle/shapeit only)* |

---

## 3. Per-Replicate Pileup (`pileup_dir/`)

One subdirectory per `{assay_type}_{rep_id}` combination.

| File | Description |
|------|-------------|
| `{assay_type}_{rep_id}/cellSNP.base.vcf.gz` | Per-replicate SNP VCF from the pileup at phased het SNP positions. |
| `{assay_type}_{rep_id}/cellSNP.samples.tsv` | Sample or barcode list for this replicate's pileup. |
| `{assay_type}_{rep_id}/cellSNP.tag.DP.mtx` | Total depth sparse matrix (SNPs × cells/samples) for this replicate. |
| `{assay_type}_{rep_id}/cellSNP.tag.AD.mtx` | Alt-allele depth sparse matrix (SNPs × cells/samples) for this replicate. |

---

## 4. Per-Assay Allele Counts (`allele_dir/{assay_type}/`)

### 4a. Feature Count Matrix (single-cell / spatial only)

| File | Description |
|------|-------------|
| `{assay_type}.h5ad` | AnnData object with a cells × genomic-bins count matrix built from gene expression (scRNA/VISIUM) or ATAC tile accessibility (scATAC), with region metadata mapped to the bin BED. |

### 4b. SNP-Level Allele Matrices (all modes)

| File | Description |
|------|-------------|
| `snps.tsv.gz` | SNP annotation table: chromosome, position, ref/alt alleles, phase, and block membership. |
| `snp.Tallele.npz` | Sparse matrix (SNPs × samples/cells): total allele count per SNP. |
| `snp.Aallele.npz` | Sparse matrix (SNPs × samples/cells): A-allele (phased haplotype 1) count per SNP. |
| `snp.Ballele.npz` | Sparse matrix (SNPs × samples/cells): B-allele (phased haplotype 2) count per SNP. |
| `sample_ids.tsv` | Sample metadata table: replicate IDs, sample types (normal/tumor). |
| `barcodes.tsv.gz` | Cell barcode or bulk sample ID list corresponding to matrix columns. |
| `unique_snp_ids.npy` | Numpy array of unique SNP identifiers used for cross-replicate alignment. *(single-cell modes only)* |
| `cell_snp_Aallele.npz` | Full cells × SNP sparse matrix of A-allele counts before binning. *(single-cell modes only)* |
| `cell_snp_Ballele.npz` | Full cells × SNP sparse matrix of B-allele counts before binning. *(single-cell modes only)* |

### 4c. Multi-SNP Block Aggregates (adaptive binning intermediate)

| File | Description |
|------|-------------|
| `multi_snp.tsv.gz` | Annotation table for multi-SNP haplotype blocks after grouping by phase and switch probabilities. |
| `multi_snp.Tallele.npz` | Total allele count matrix aggregated to multi-SNP block level. |
| `multi_snp.Aallele.npz` | A-allele count matrix at multi-SNP block level. |
| `multi_snp.Ballele.npz` | B-allele count matrix at multi-SNP block level. |

### 4d. Final Genomic Bins ("bb" = BAF blocks) — primary outputs

| File | Description |
|------|-------------|
| `bb.tsv.gz` | Genomic bin annotation table: chromosome, start/end, SNP count, gene/tile annotations. One row per bin. |
| `bb.Tallele.npz` | Sparse matrix (bins × samples/cells): total allele count per bin. |
| `bb.Aallele.npz` | Sparse matrix (bins × samples/cells): A-allele count per bin. |
| `bb.Ballele.npz` | Sparse matrix (bins × samples/cells): B-allele count per bin. |
| `bb.bed.gz` | BED file of genomic bins; used as input regions for `mosdepth` coverage computation. |
| `bb.depth.npz` | Read depth matrix (bins × samples) from `mosdepth`. *(bulk assays only)* |
| `bb.rdr.npz` | Read depth ratio (RDR) matrix normalized to a diploid baseline. *(bulk assays only)* |
| `corr_factors.tsv.gz` | GC-content (and optionally mappability) correction factors per bin used to compute RDR. *(bulk assays only)* |
| `out_mosdepth/{rep_id}.regions.bed.gz` | Raw per-replicate `mosdepth` output: read depth over each bin. *(bulk assays only)* |

### 4e. CNV Segmentation Outputs (`copytyping_preprocess` mode only)

| File | Description |
|------|-------------|
| `cnv_segments.tsv` | CNV segment table annotated with integer total and allele-specific copy numbers from an input segmentation. |
| `X_count.npz` | Sparse matrix (cells × segments): total allele count per cell per CNV segment. |
| `Y_count.npz` | Sparse matrix (cells × segments): A-allele count per cell per CNV segment. |
| `D_count.npz` | Sparse matrix (cells × segments): read depth count per cell per CNV segment. |

---

## 5. Column Reference for TSV Files

### `snps.tsv.gz`

| Column | Description |
|--------|-------------|
| `#CHR` | Chromosome name. |
| `POS` | 1-based SNP position (right end, i.e. `POS0 + 1`). |
| `POS0` | 0-based SNP position. |
| `START` | Start of the SNP's assigned region interval (from `region_bed`). |
| `END` | End of the SNP's assigned region interval. |
| `GT` | Genotype string (e.g. `0|1`, `1|0`, `1|1`). |
| `PHASE` | Phased allele orientation: `1` if A-allele is ALT, `-1` if A-allele is REF. |
| `region_id` | Integer ID of the genomic region (from `region_bed`) this SNP falls in. |
| `feature_id` | Integer ID of the gene/tile feature this SNP is assigned to (from `gtf_file`). |

### `multi_snp.tsv.gz`

| Column | Description |
|--------|-------------|
| `#CHR` | Chromosome name. |
| `START` | Genomic start of the bin (min over member SNP region starts). |
| `END` | Genomic end of the bin (max over member SNP region ends). |
| `START0` | 0-based position of the leftmost SNP in the bin. |
| `END0` | 1-based position of the rightmost SNP in the bin. |
| `switchprobs` | Haplotype-switch probability at the first SNP of this bin (inter-bin boundary). Present only when a genetic map or PS-based switch prob is available. |
| `region_id` | Region ID shared by all SNPs in this bin. |
| `#SNPS` | Number of SNPs aggregated into this bin. |
| `BLOCKSIZE` | Genomic span of the bin in bp (`END - START`). |
| `multi_id` | Integer bin index (0-based, global across all chromosomes). |

### `bb.tsv.gz`

Same structure as `multi_snp.tsv.gz` but with additional grouping columns and a different bin index:

| Column | Description |
|--------|-------------|
| `#CHR` | Chromosome name. |
| `START` | Genomic start of the bin. |
| `END` | Genomic end of the bin. |
| `START0` | 0-based position of the leftmost SNP. |
| `END0` | 1-based position of the rightmost SNP. |
| `switchprobs` | Inter-bin haplotype-switch probability at the bin boundary. |
| `region_id` | Region ID shared by all SNPs in this bin. |
| `PS` | Phase-set ID; SNPs in the same bin share the same haplotype phase. |
| `binom_id` | Binomial-test boundary ID; breaks bins at significant AF shifts. Present only when `binom_test: true`. |
| `#SNPS` | Number of SNPs in this bin. |
| `BLOCKSIZE` | Genomic span in bp. |
| `bb_id` | Integer bin index (0-based, global across all chromosomes). |

### `sample_ids.tsv`

| Column | Description |
|--------|-------------|
| `SAMPLE` | Full sample label (`{SAMPLE_NAME}_{REP_ID}`). |
| `SAMPLE_NAME` | Patient-level sample name (same for all rows). |
| `REP_ID` | Replicate identifier. |
| `sample_type` | `normal` or `tumor`. |

### `barcodes.tsv.gz`

Single column, no header. Each row is one cell barcode (or bulk sample label) in the form `{BARCODE}_{REP_ID}`. Row order matches the columns of the allele count matrices.

### `corr_factors.tsv.gz`

Contains all columns from `bb.tsv.gz` plus:

| Column | Description |
|--------|-------------|
| `GC` | GC-content fraction of the bin computed from the reference genome. |
| `MAP` | Mean mappability score of the bin (`1.0` if no mappability file is provided). |

### `cnv_segments.tsv`

Inherits all columns from the input HATCHet `.ucn` file (typically `#CHR`, `START`, `END`, per-clone `cn_*` and `u_*` columns), plus:

| Column | Description |
|--------|-------------|
| `seg_id` | Integer segment index (0-based). |
| `CNP` | Semicolon-separated total copy numbers across all clones. |
| `PROPS` | Semicolon-separated clone proportions (usage values) across all clones. |
