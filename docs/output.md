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
