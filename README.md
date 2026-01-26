## Unified Single-cell/Visium Genotyping+Phasing+Pileup Workflow
```
mamba env create -f ./environment.yaml -p /path/to/envs/genotyping_env
conda activate /path/to/envs/genotyping_env
snakemake -p --cores 1 -s /path/to/workflow/Snakefile --configfile config/config.yaml --directory <output>
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
3. dbSNPv155
    1. https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/liftover/*

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

Some requirements on sample files:
* All rows should have same `SAMPLE` field, replicates can be distinguished based on replicate ID `REP_ID`
* Paired scMultiome can be described by two rows with same `REP_ID` and two `DATA_TYPE` values: `scRNA` and `scATAC`.
* For mode 1, specify `bulkDNA` under `DATA_TYPE`, and use `REP_ID`=`normal` to label the matched-normal sample if available.
* For mode 2, `DATA_TYPE` can be chosen from `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime`.

### Output


### TODO
1. Detecting clonal LOH, from Numbat
```
In samples with high tumor purity (e.g., tumor cell lines) without matched normal cells, heterozygous SNPs are challenging to identify in regions of LoH, leading to decreased power of CNV detection. Regions of clonal LoH have decreased SNP density and can be identified by a specialized HMM (see detect_clonal_loh in function reference).
```
2. mapping resolutions, 2um to Xum?
3. divides post het SNPs into haplotype blocks, compute RDR, then aggregate to form phased blocks.
