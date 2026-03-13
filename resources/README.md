# Resources

External data required by the pipeline.

---

## SNP Panels

VCF format. Set via `snp_panel` or `snp_targets` in config.

| Panel | hg38 |
|-------|------|
| 1kGP phase3 AF>=5e-2 (~92 MB) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) |
| 1kGP phase3 AF>=5e-4 (~568 MB) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz) |
| 1kGP n=3,202 (high-coverage) | [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) — use [`scripts/process_1kGP_3202_panel.sh`](scripts/process_1kGP_3202_panel.sh) to prepare |

### Building `snp_targets` from any panel VCF

The `genotype_snps_bulk` rule requires `config["snp_targets"]` — a directory of per-chromosome position files. To build these from any SNP panel VCF:

```bash
bash resources/scripts/build_snp_targets.sh /path/to/snp_panel.vcf.gz /path/to/snp_targets
```

This produces `target.chr{1..22,X}.pos.gz` + `.tbi` index files. Existing chromosomes are skipped (idempotent). Requires `bcftools`, `bgzip`, `tabix`.

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
