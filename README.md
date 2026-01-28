## Unified Single-cell/Visium Genotyping+Phasing+Pileup Workflow
```
mamba env create -f ./environment.yaml -p /path/to/envs/genotyping_env
conda activate /path/to/envs/genotyping_env
snakemake --cores <num_cores> \
    -s /path/to/workflow/Snakefile \
    --report report.html \
    --configfile config/config.yaml \
    --directory <output>
```

### Dependencies
1. samtools
2. bcftools
3. cellsnp-lite
4. HetDetect `https://github.com/raphael-group/hetdetect.git`
5. Eagle2 `https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz`
6. Shapeit5 `https://github.com/odelaneau/shapeit5.git`, `https://github.com/odelaneau/shapeit5/releases`

### Aux files
#### SNP panel files
1. https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
    * used by Numbat, ~92MB
2. https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz
    * used by CalicoST, ~568MB
3. dbSNPv157 (for WGS data)
    * hg19: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
    * hg38: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz and https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

##### rename dbSNP SNP panel chromosomes
NCBI dbSNP uses refseq version, rename it to UCSC version.
```sh
bcftools view -r NC_000001.11,NC_000002.12,NC_000003.12,NC_000004.12,NC_000005.10,NC_000006.12,NC_000007.14,NC_000008.11,NC_000009.12,NC_000010.11,NC_000011.10,NC_000012.12,NC_000013.11,NC_000014.9,NC_000015.10,NC_000016.10,NC_000017.11,NC_000018.10,NC_000019.10,NC_000020.11,NC_000021.9,NC_000022.11,NC_000023.11,NC_000024.10,NC_012920.1 -Ou GCF_000001405.40.gz \
    | bcftools annotate --rename-chrs /path/to/rename_chrs.refseq2ucsc.tsv -Ou \
    | bcftools view -v snps -Oz -o dbsnp157.hg38.biallelic.snps.vcf.gz

bcftools index dbsnp157.hg38.biallelic.snps.vcf.gz
```

#### phasing panel files
1. 1000GP HG38 phasing panel: http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip
2. gnomAD HGDP + 1KG panel (n=4,099): download the reference files using gsutil: gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes
3. TOPMed imputation server (n=97,256)

### Running options
User can either run with mode 1 (bulk WGS/WES genotyping) or mode 2 (single-cell/visium genotyping). Refer to mode1 and mode2 DAG workflow for all involved rules.
* For mode 1, sample file should only contain bulk samples.
* For mode 2, sample file only only contain non-bulk samples.
* In either mode, if reference het SNP file genotyped from same patient is provided, genotyping and phasign step are skipped, only pile-up and postprocessing will run.

#### Sample file
Sample file is tsv-format and has following columns to describe all samples from same patient:
* (Required columns): SAMPLE, REP_ID, DATA_TYPE, PATH_to_bam
* (Optional columns): PATH_to_barcodes, PATH_to_10x_ranger

Requirements on sample files:
* All rows should have same `SAMPLE` field, replicates can be distinguished based on replicate ID `REP_ID`
* Paired scMultiome can be described by two rows with same `REP_ID` and two `DATA_TYPE` values: `scRNA` and `scATAC`.
* For mode 1, specify `bulkDNA` under `DATA_TYPE`, and use `REP_ID`=`normal` to label the matched-normal sample if available.
* For mode 2, `DATA_TYPE` can be chosen from `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime`.
* provide two separate sample files for mode 1 and mode 2, respectively. Don't interleave bulk and non-bulk samples in same sample file.
* when both bulk and non-bulk samples from same patient are provided, user should first run mode 1 using bulk samples, set `ref_snp_file=phase/phased_snps.vcf.gz` in config file, and run mode 2 on non-bulk samples.

### Output
* `snps/<chromosome>.vcf.gz*`: annotated bi-allelic Het/Hom-Alt SNPs.
* `phase/phased_snps.vcf.gz*`: phased Het SNPs.
* `allele/<data_type>(_<rep_id>)?/sample_ids.tsv`: sample IDs.
* `allele/<data_type>(_<rep_id>)?/snp_info.tsv.gz`: Per-SNP positional and genotype information.
* `allele/<data_type>(_<rep_id>)?/snp_matrix.dp.npz`: SNP by samples read-depth matrix, generated if mode 1.
* `allele/<data_type>(_<rep_id>)?/snp_matrix.tot.npz`: SNP by samples/barcodes total-allele matrix.
* `allele/<data_type>(_<rep_id>)?/snp_matrix.ref.npz`: SNP by samples/barcodes ref-allele matrix.
* `allele/<data_type>(_<rep_id>)?/snp_matrix.alt.npz`: SNP by samples/barcodes alt-allele matrix.
* `allele/<data_type>(_<rep_id>)?/snp_matrix.A.npz`: SNP by samples/barcodes A-allele matrix, generated if phased.
* `allele/<data_type>(_<rep_id>)?/snp_matrix.B.npz`: SNP by samples/barcodes B-allele matrix, generated if phased.

#### Optional output
* `phase/genetic_map.tsv.gz`: parsed genetic map file, required by HATCHet3.

#### Legacy output
* `allele/<data_type>(_<rep_id>)?/unique_snp_ids.npy`: per-SNP `<chromosome>_<position>` array.
* `allele/<data_type>(_<rep_id>)?/cell_snp_Aallele.npz`: SNP by samples/barcodes A-allele matrix, generated if phased.
* `allele/<data_type>(_<rep_id>)?/cell_snp_Ballele.npz`: SNP by samples/barcodes B-allele matrix, generated if phased.

### TODO
1. Detecting clonal LOH, from Numbat
```
In samples with high tumor purity (e.g., tumor cell lines) without matched normal cells, heterozygous SNPs are challenging to identify in regions of LoH, leading to decreased power of CNV detection. Regions of clonal LoH have decreased SNP density and can be identified by a specialized HMM (see detect_clonal_loh in function reference).
```
2. mapping resolutions, 2um to Xum?
3. divides post het SNPs into haplotype blocks, compute RDR, then aggregate to form phased blocks.
4. we assumed all BAM files have chr-prefix, and we preserve chr-prefix in all processed files.
