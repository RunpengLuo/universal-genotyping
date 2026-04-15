# Resources

External data required by the pipeline.

---

## SNP Panels

VCF format. Set via `snp_panel` or `snp_targets` in config.

| Panel | Download |
|-------|----------|
| 1kGP phase3 AF>=5e-2 (~92 MB, hg38) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) |
| 1kGP phase3 AF>=5e-4 (~568 MB, hg38) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz) |
| 1kGP n=3,202 (hg38) | [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) — use [`scripts/process_1kGP_3202_panel.sh --ref hg38`](scripts/process_1kGP_3202_panel.sh) |
| 1kGP n=3,202 (chm13v2.0, biallelic) | [S3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/) — use [`scripts/process_1kGP_3202_panel.sh --ref chm13v2`](scripts/process_1kGP_3202_panel.sh) |

### Building `snp_targets` from any panel VCF

The `genotype_snps_bulk` rule requires `config["snp_targets"]` — a directory of per-chromosome position files. To build these from any SNP panel VCF:

```bash
bash resources/scripts/build_snp_targets.sh /path/to/snp_panel.vcf.gz /path/to/snp_targets
```

This produces `target.chr{1..22,X}.pos.gz` + `.tbi` index files. Existing chromosomes are skipped (idempotent). Requires `bcftools`, `bgzip`, `tabix`.

---

## Phasing Panels

BCF format, one per chromosome. Set via `phasing_panel` in config.

| Panel | Download |
|-------|----------|
| 1kGP phase3 (n=2,504, hg38) | [download](http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip) |
| 1kGP phase3 (n=3,202, hg38) | produced by `process_1kGP_3202_panel.sh --ref hg38` (see SNP Panels) |
| 1kGP n=3,202 (chm13v2.0) | produced by `process_1kGP_3202_panel.sh --ref chm13v2` (see SNP Panels). Phased with SHAPEIT5 v1.1 + T2T-native maps ([phasing_T2T](https://github.com/JosephLalli/phasing_T2T)) |
| gnomAD HGDP+1KG (n=4,099, hg38) | `gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes` |
| TOPMed (n=97,256, hg38) | via [imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov) |

### Genetic Maps

- **hg38:** bundled with Eagle2 (`tables/`) and SHAPEIT5 (`resources/maps/`). Set `phaser_dir` to the phasing program root directory.
- **chm13v2:** download [T2T-native scaled maps](https://github.com/JosephLalli/phasing_T2T/tree/main/resources/recombination_maps/t2t_native_scaled_maps) and place under `{phaser_dir}`:
  - **SHAPEIT5:** copy maps to `{phaser_dir}/resources/maps/chm13v2/` (files named `chr{N}.t2t.scaled.gmap.gz`)
  - **Eagle2:** convert with [`scripts/convert_gmap_to_eagle.py`](scripts/convert_gmap_to_eagle.py) and place output at `{phaser_dir}/tables/genetic_map_chm13v2_withX.txt.gz`

---

## Gene Annotation (GTF)

Set via `gtf_file` in config.

| Source | Download |
|--------|----------|
| GENCODE v38 (hg38) | `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz` |
| 10x GRCh38-2024-A (hg38) | `curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` → `genes/genes.gtf.gz` |
| UCSC ncbiRefSeq (chm13v2.0, chr-style) | `wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz` |
| NCBI RefSeq (chm13v2.0, accession-style) | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/) — `GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz` |

Build 10x Cell Ranger ARC reference with [`scripts/build_cellranger_arc_ref_chm13v2.sh`](scripts/build_cellranger_arc_ref_chm13v2.sh) ([10x guide](https://kb.10xgenomics.com/s/article/29207065679501-Building-a-Custom-T2T-reference-for-Cell-Ranger-ARC)).

---

## Window BED (Bias Correction)

Set via `window_bed` in config. A pre-built file is at `data/windows.1kbp.hg38.bed.gz` with columns for GC content, mappability, and replication timing per 1 kb window.

To build from scratch:
- **WGS:** [`scripts/build_wgs_window_bed.py`](scripts/build_wgs_window_bed.py) — fixed-size tiling with region/blacklist filtering.
- **WES:** [`scripts/build_wes_window_bed.py`](scripts/build_wes_window_bed.py) — adaptive tiling of WES capture targets via `--wes_targets_bed`.

### Mappability track

Optional, used by `build_wgs_window_bed.py --mappability_bed` to add the `MAP` column.

- [k100.Umap.MultiTrackMappability.bw](http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw) — bigWig format, convert to BED with `bigWigToBedGraph` (UCSC tools).

### Replication timing (Repli-seq)

Optional, used by `build_wgs_window_bed.py --repliseq` to add the `REPLI` column. The script automatically downloads 16 ENCODE Repli-seq WaveSignal bigWig files from [UCSC](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/) (hg19), converts via `bigWigToBedGraph`, and lifts to hg38 using [hg19ToHg38.over.chain.gz](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz). Requires `bigWigToBedGraph` and `liftOver` (UCSC tools).

### WES exon capture targets

For WES mode, `build_wes_window_bed.py --wes_targets_bed` requires a vendor exon capture BED. Example (IDT xGen):

- [xgen-exome-research-panel-targets-hg38.bb](https://hgdownload.soe.ucsc.edu/gbdb/hg38/exomeProbesets/xgen-exome-research-panel-targets-hg38.bb) — bigBed format, convert to BED with `bigBedToBed` (UCSC tools).

---

## Blacklist BED

Set via `blacklist_bed` in config. Pre-built: `data/hg38-blacklist.v2.bed.gz` ([ENCODE blacklist v2](https://github.com/Boyle-Lab/Blacklist)).

---

## Phasing Tools

One required for `bulk_genotyping` and `single_cell_genotyping` modes.

| Tool | Source |
|------|--------|
| Eagle2 | [download](https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz) |
| SHAPEIT5 | [GitHub](https://github.com/odelaneau/shapeit5) |
| LongPhase | [GitHub](https://github.com/twolinin/longphase) |
