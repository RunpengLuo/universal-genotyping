# Resources

External data files required to run the pipeline.

---

## Reference Genome

Required by alignment, genotyping, and blacklist construction (`reference` in config).

| Version | Contigs | Use case | Download |
|---------|---------|----------|----------|
| 10x Genomics GRCh38 | chr1–22, X, Y only | Single-cell pipelines (Cell Ranger, spaceranger) | `curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` → `fasta/genome.fa` |
| NCBI GRCh38 analysis set | chr1–22, X, Y + decoys, HLA, EBV (no ALT contigs) | Bulk WGS/WES alignment, mappability generation | see below |

```bash
# NCBI GRCh38 analysis set (no ALT contigs, with decoys + HLA)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

The pipeline only uses chr1–22, X, Y for SNP calling and phasing. Extra contigs in the NCBI set act as decoy targets during alignment to reduce mismapping.

---

## SNP Panels

Used for genotyping (`snp_panel` or `snp_targets` in config). VCF format, one file per genome or per chromosome.

| Panel | hg19 | hg38 | Notes |
|-------|------|------|-------|
| 1kGP phase3 AF≥5e-2 | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) | ~92 MB; used by Numbat |
| 1kGP phase3 AF≥5e-4 | — | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz) | ~568 MB; used by CalicoST |
| 1kGP phase3 n=3,202 (high-coverage) | — | [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) · [README](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf) | Per-chromosome VCFs; also doubles as phasing panel. Use [`resources/scripts/process_1kGP_3202_panel.sh`](../resources/scripts/process_1kGP_3202_panel.sh) to download and prepare all outputs. |

---

## Phasing Panels

Used for population-based phasing (`phasing_panel` in config). BCF format, one file per chromosome.

| Panel | hg19 | hg38 | Notes |
|-------|------|------|-------|
| 1kGP phase3 (n=2,504) | [download](http://pklab.med.harvard.edu/teng/data/1000G_hg19.zip) | [download](http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip) | Commonly used with Eagle/SHAPEIT |
| 1kGP phase3 (n=3,202, with trios) | — | see SNP Panels above | Higher coverage; includes 698 trios. BCF outputs produced by [`resources/scripts/process_1kGP_3202_panel.sh`](../resources/scripts/process_1kGP_3202_panel.sh). |
| gnomAD HGDP + 1KG (n=4,099) | — | `gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes` (gsutil) | Diverse ancestry panel |
| TOPMed (n=97,256) | — | via [imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov) | Largest panel; requires server access |

---

## Gene Annotation (GTF)

Used for feature assignment (`gtf_file` in config).

| Source | Download |
|--------|----------|
| GENCODE v38 (hg38) | `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz` |
| 10x Genomics GRCh38-2024-A | `curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` → `genes/genes.gtf.gz` |

---

## Utilities

| Tool | Download |
|------|----------|
| subset-bam (10x) | [download](https://github.com/10XGenomics/subset-bam/releases/download/v1.1.0/subset-bam_linux) |

---

## Blacklist BED

Used for filtering unreliable genomic regions (`blacklist_file` in config).

### Building the blacklist

The script [`resources/scripts/build_blacklist.hg38.py`](../resources/scripts/build_blacklist.hg38.py) generates a merged blacklist BED by:

1. Running hmmcopy_utils (`generateMap.pl` + `mapCounter`) to compute per-base mappability from a reference FASTA
2. Filtering bins with low average mappability
3. Downloading UCSC repeat/segdup tables (genomicSuperDups, simpleRepeat)
4. Merging all intervals with pyranges into a single gzipped BED

**Dependencies:** [hmmcopy_utils](https://github.com/shahcompbio/hmmcopy_utils) (`generateMap.pl`, `mapCounter`) and [bowtie](https://bowtie-bio.sourceforge.net/index.shtml).

**Quick start:**

```bash
python resources/scripts/build_blacklist.hg38.py \
  --reference /path/to/hg38.fa \
  --out_file /path/to/blacklist.hg38.bed.gz \
  --build_index --threads 8

# Keep intermediate files (mappability BigWig, etc.) for inspection
python resources/scripts/build_blacklist.hg38.py \
  --reference /path/to/hg38.fa \
  --out_file /path/to/blacklist.hg38.bed.gz \
  --build_index --threads 8 \
  --work_dir /path/to/intermediates
```

**Choosing `--read_length`:** Use the default **150** for both standard short-read WGS/WES and PacBio HiFi workflows. `generateMap.pl` aligns synthetic k-mers with bowtie (a short-read aligner), so passing long-read lengths (10 kb+) is not supported. The blacklist also includes UCSC segdups and simple repeats regardless of read length.

Run `python resources/scripts/build_blacklist.hg38.py --help` for full usage.

---

## Phasing Tools

One phasing tool is required for `bulk_genotyping` and `single_cell_genotyping` modes. Python/conda dependencies are in `environment.yaml`.

| Tool | Source |
|------|--------|
| Eagle2 | [download](https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz) |
| SHAPEIT5 | [GitHub](https://github.com/odelaneau/shapeit5) |
| LongPhase | [GitHub](https://github.com/twolinin/longphase) |

