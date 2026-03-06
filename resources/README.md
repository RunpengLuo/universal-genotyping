# Resources

External data required by the pipeline.

---

## Reference Genome

| Version | Use case | Download |
|---------|----------|----------|
| 10x GRCh38 | Single-cell (Cell Ranger) | `curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` → `fasta/genome.fa` |
| NCBI GRCh38 analysis set | Bulk WGS/WES | `wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz` |

The pipeline only uses chr1-22, X, Y for SNP calling and phasing.

---

## SNP Panels

VCF format. Set via `snp_panel` or `snp_targets` in config.

| Panel | hg38 |
|-------|------|
| 1kGP phase3 AF>=5e-2 (~92 MB) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) |
| 1kGP phase3 AF>=5e-4 (~568 MB) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz) |
| 1kGP n=3,202 (high-coverage) | [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) — use [`scripts/process_1kGP_3202_panel.sh`](scripts/process_1kGP_3202_panel.sh) to prepare |

---

## Phasing Panels

BCF format, one per chromosome. Set via `phasing_panel` in config.

| Panel | hg38 |
|-------|------|
| 1kGP phase3 (n=2,504) | [download](http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip) |
| 1kGP phase3 (n=3,202) | see SNP Panels above |
| gnomAD HGDP+1KG (n=4,099) | `gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes` |
| TOPMed (n=97,256) | via [imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov) |

---

## Gene Annotation (GTF)

Set via `gtf_file` in config.

| Source | Download |
|--------|----------|
| GENCODE v38 | `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz` |
| 10x GRCh38-2024-A | `curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` → `genes/genes.gtf.gz` |

---

## Bias Correction BED

Set via `params_compute_rdr.bias_bed` in config. A pre-built file is at `data/gc_map_repli.1kbp.hg38.bed.gz` with columns for GC content, mappability, and replication timing per 1 kb window.

To build from scratch, see `scripts/build_window_bed.py`.

---

## Blacklist BED

Set via `blacklist_bed` in config. Pre-built: `data/hg38-blacklist.v2.bed.gz` ([ENCODE blacklist v2](https://github.com/Boyle-Lab/Blacklist)).

https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed

https://hgdownload.soe.ucsc.edu/gbdb/hg38/exomeProbesets/xgen-exome-research-panel-targets-hg38.bb
---

## Phasing Tools

One required for `bulk_genotyping` and `single_cell_genotyping` modes.

| Tool | Source |
|------|--------|
| Eagle2 | [download](https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz) |
| SHAPEIT5 | [GitHub](https://github.com/odelaneau/shapeit5) |
| LongPhase | [GitHub](https://github.com/twolinin/longphase) |
