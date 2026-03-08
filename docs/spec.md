# Input / Output Reference

## Sample File

TSV with one row per replicate. Template at `config/samples.tsv`.

| Column | Required | Description |
|--------|----------|-------------|
| `SAMPLE` | Yes | Patient ID (matched by `sample_id` in config). |
| `REP_ID` | Yes | Unique per row; multiome pairs share a `REP_ID`. |
| `assay_type` | Yes | `bulkWGS`, `bulkWES`, `scATAC`, `scRNA`, `VISIUM`, or `VISIUM3prime`. |
| `sample_type` | Yes | `normal` or `tumor`. |
| `PATH_to_bam` | Yes | Path to `.bam` file. |
| `PATH_to_barcodes` | Non-bulk | Path to `barcodes.tsv.gz`. |
| `PATH_to_10x_ranger` | Non-bulk | Path to Cell Ranger / Space Ranger `outs/` directory. |
| `PATH_to_ref_annotations` | Optional | Annotation TSV with `BARCODE` and `<ref_label>` columns. |

---

## Outputs

All output directories (`snp_dir`, `phase_dir`, `pileup_dir`, `allele_dir`, `bb_dir`, `qc_dir`, `log_dir`) are set in `config.yaml`.

### SNP Genotyping (`snp_dir/`)

- `chr{chrname}.vcf.gz` — per-chromosome VCF of bi-allelic SNPs.
- `pseudobulk_{modality}/cellSNP.*` — pseudobulk pileup (single-cell only).

### Phasing (`phase_dir/`)

- `chr{chrname}.vcf.gz` — per-chromosome phased VCF.
- `phased_het_snps.vcf.gz(.tbi)` — genome-wide phased het SNPs.
- `genetic_map.tsv.gz` — cM positions per SNP (eagle/shapeit only).

### Pileup (`pileup_dir/`)

One subdirectory per `{assay_type}_{rep_id}` with cellsnp-lite output.

### Allele Counts (`allele_dir/{assay_type}/`)

- `snps.tsv.gz` — SNP annotations (chr, pos, ref/alt, phase, block).
- `snp.{Tallele,Aallele,Ballele}.npz` — sparse allele count matrices (SNPs x samples/cells).
- `barcodes.tsv.gz` — cell barcode list (single-cell).
- `{assay_type}.h5ad` — AnnData with cells x bins (single-cell).

### Final Bins (`bb_dir/{assay_type}/`)

**Bulk (`bulk_genotyping`):**
- `bb.tsv.gz` — bin annotations.
- `bb.{Tallele,Aallele,Ballele,baf,depth,rdr}.npz` — allele, depth, and RDR matrices.

**Single-cell (`single_cell_genotyping`):**
- `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele}.npz`

**Copytyping (`copytyping_preprocess`):**
- `cnv_segments.tsv`, `{X,Y,D}_count.npz`

### QC (`qc_dir/`)

Plots from bias correction and binning steps, organized under `qc_dir/{assay_type}/`.

---

## Key TSV Columns

### `snps.tsv.gz`

`#CHR`, `POS`, `POS0`, `START`, `END`, `GT`, `PHASE` (0 = B-allele is ALT, 1 = B-allele is REF), `region_id`, `feature_id`, `feature_type` (exon/intron/intergenic).

### `bb.tsv.gz`

`#CHR`, `START`, `END`, `START0`, `END0`, `switchprobs`, `region_id`, `PS`, `#SNPS`, `BLOCKSIZE`, `feature_id`, `feature_type`, `bb_id`.

### `sample_ids.tsv`

`SAMPLE`, `SAMPLE_NAME`, `REP_ID`, `sample_type`.

### `barcodes.tsv.gz`

Single column, no header. Each row: `{BARCODE}_{REP_ID}`.
