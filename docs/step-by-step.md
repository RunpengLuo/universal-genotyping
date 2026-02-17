# Step-by-Step Guide

All three modes share the same command structure — only `workflow_mode` in `config/config.yaml` and the sample file differ. Fill in the paths in `config/config.yaml` (reference genome, SNP panel, phasing panel, tool binaries, output directories, etc.) before running any mode.

```sh
snakemake --cores <N> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output_dir> \
    --config sample_file=/path/to/samples.tsv
```

---

## Mode 1: `bulk_genotyping`

Use this mode when you have bulk WGS or WES BAM files. The pipeline genotypes bi-allelic SNPs from the matched-normal BAM via `bcftools`, phases them using the configured phaser (Eagle, SHAPEIT, or LongPhase), pileups allele counts per sample, and computes GC-corrected read-depth ratios (RDR).

1. Prepare a sample sheet listing only bulk samples (`assay_type` ∈ `{bulkWGS, bulkWES}`). Include both normal and tumor rows, setting `sample_type` accordingly.
2. Set `workflow_mode: bulk_genotyping` in `config/config.yaml`.
3. Run the pipeline. Primary outputs in `allele_dir/{assay_type}/`: `bb.tsv.gz`, `bb.Tallele.npz`, `bb.Aallele.npz`, `bb.Ballele.npz`, `bb.depth.npz`, `bb.rdr.npz`.

---

## Mode 2: `single_cell_genotyping`

Use this mode when you have single-cell or spatial data (scRNA, scATAC, VISIUM, VISIUM3prime). SNPs are genotyped from a pseudobulk pileup across all BAMs using `cellsnp-lite`, phased, then piled up per cell. Gene expression or ATAC accessibility is aggregated into the same genomic bins.

1. Prepare a sample sheet listing only non-bulk samples (`assay_type` ∈ `{scRNA, scATAC, VISIUM, VISIUM3prime}`). For paired scMultiome, give the scRNA and scATAC rows the same `REP_ID`. Fill in `PATH_to_barcodes` and `PATH_to_10x_ranger` for every row.
2. Set `workflow_mode: single_cell_genotyping` in `config/config.yaml`.
3. Run the pipeline. Primary outputs in `allele_dir/{assay_type}/`: `bb.tsv.gz`, `bb.Tallele.npz`, `bb.Aallele.npz`, `bb.Ballele.npz`, `barcodes.tsv.gz`, `{assay_type}.h5ad`.

**Tip:** If you already ran Mode 1 on matched bulk data from the same patient, you can skip re-genotyping by setting `het_snp_vcf: <output_dir>/phase/phased_het_snps.vcf.gz` in the config. Only pileup and postprocessing will run.

---

## Mode 3: `copytyping_preprocess`

Use this mode when phased het SNPs and a CNV profile (e.g., from HATCHet3) are already available. The pipeline skips genotyping and phasing entirely, and instead aggregates per-cell allele counts directly onto the provided CNV segments.

1. Prepare a sample sheet listing non-bulk samples (same as Mode 2).
2. Set `workflow_mode: copytyping_preprocess` in `config/config.yaml`, and provide:
   - `het_snp_vcf`: path to phased het SNP VCF (e.g., from a prior Mode 1 run).
   - `seg_ucn`: path to HATCHet3 segment-level copy-number file.
3. Run the pipeline. Primary outputs in `allele_dir/{assay_type}/`: `cnv_segments.tsv`, `X_count.npz`, `Y_count.npz`, `D_count.npz`.
