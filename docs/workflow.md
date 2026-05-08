# Workflow

Three workflow modes, set via `config["workflow_mode"]`. Each mode runs a subset of the pipeline stages.

---

## `bulk_genotyping`

Assays: `bulkWGS`, `bulkWGS-lr`, `bulkWES`

| Step | Rule | Script / Tool |
|------|------|---------------|
| 1. Genotype SNPs | `genotype_snps_bulk` | bcftools mpileup + call |
| 2. Phase SNPs (per chr) | `phase_snps_{eagle,shapeit,longphase}` | eagle2 / shapeit5 / longphase |
| 3. Concat phased VCFs | `concat_and_extract_phased_het_snps` | bcftools concat + view |
| 4. Parse genetic map | `parse_genetic_map` | `scripts/parse_genetic_map.py` |
| 5. Pileup at het SNPs | `pileup_snps_bulk_mode1b` | cellsnp-lite |
| 6. Phase and concat | `phase_and_concat_bulk` | `scripts/phase_and_concat.py` |
| 7. Compute read depth | `run_mosdepth` | mosdepth |
| 8. Bias correction | `rd_correct` | `scripts/rd_correct.py` |
| 9. Adaptive binning | `combine_counts` | `scripts/combine_counts.py` |

**Outputs** (`bb_dir`): `bb.{tsv.gz,Tallele,Aallele,Ballele,baf,rdr,depth}.npz`, `sample_ids.tsv`

---

## `single_cell_genotyping`

Assays: `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime`

| Step | Rule | Script / Tool |
|------|------|---------------|
| 1. Pseudobulk genotyping | `genotype_snps_pseudobulk_mode1b` | cellsnp-lite |
| 2. Annotate SNPs | `annotate_snps_pseudobulk` | `scripts/annotate_snps_pseudobulk.py` |
| 3. Phase SNPs (per chr) | `phase_snps_{eagle,shapeit,longphase}` | eagle2 / shapeit5 / longphase |
| 4. Concat phased VCFs | `concat_and_extract_phased_het_snps` | bcftools concat + view |
| 5. Parse genetic map | `parse_genetic_map` | `scripts/parse_genetic_map.py` |
| 6. Single-cell pileup | `pileup_snps_single_cell_mode1a` | cellsnp-lite |
| 7. Build AnnData | `process_rna_anndata` / `process_atac_fragments` | `scripts/process_rna_anndata.py` / `scripts/process_atac_fragments.py` |
| 8. Phase and concat | `phase_and_concat_single_cell` | `scripts/phase_and_concat.py` |
| 9. Adaptive binning | `combine_counts_nonbulk` | `scripts/combine_counts_nonbulk.py` |

**Outputs** (`bb_dir`): `bb.{tsv.gz,Tallele,Aallele,Ballele}.npz`, `sample_ids.tsv`

---

## `copytyping_preprocess`

Assays: `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime` (requires pre-computed het SNP VCF)

| Step | Rule | Script / Tool |
|------|------|---------------|
| 1. Single-cell pileup | `pileup_snps_single_cell_mode1a` | cellsnp-lite |
| 2. Build AnnData | `process_rna_anndata` / `process_atac_fragments` | `scripts/process_rna_anndata.py` / `scripts/process_atac_fragments.py` |
| 3. Phase and concat | `phase_and_concat_single_cell` | `scripts/phase_and_concat.py` |
| 4. Adaptive binning | `combine_counts_nonbulk` | `scripts/combine_counts_nonbulk.py` |
| 5. CNV segmentation | `cnv_segmentation` | `scripts/cnv_segmentation.py` |

**Outputs** (`bb_dir`): `cnv_segments.tsv`, `bb.{Tallele,Aallele,Ballele,Xcount}.npz`, `barcodes.tsv.gz`, `sample_ids.tsv`
