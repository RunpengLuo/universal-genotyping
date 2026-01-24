## Unified Single-cell/Visium Genotyping+Phasing+Pileup Workflow
```
mamba env create -f ./environment.yaml -p /path/to/envs/genotyping_env
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
1. http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip

gnomAD HGDP + 1KG panel (n=4,099). You can download the reference files using gsutil: gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes.

TOPMed panel (n=97,256). You can upload your VCF to the TOPMed imputation server.

### Mode 1. bulk sample
If reference Het SNP file from same patient is provided, skip step 1-2.
1. genotype SNPs via bcftools mpileup+call and keep bi-allelic SNPs.
    1. normal sample or tumor sample.
2. phase SNPs via population-based phasing methods
3. pile-up allele counts via cellsnp-lite mode 1b.
4. postprocess allele count matrices.

### Mode 2. single-cell/Visium sample
If reference Het SNP file from same patient is provided, skip step 1-2.
1. genotype SNPs via cellsnp-lite on pseudobulk sample.
2. phase SNPs via population-based phasing methods
3. pile-up allele counts via cellsnp-lite mode 1a.
4. postprocess allele count matrices.

### Sample file
Sample file has following modes:
1. multiome (scRNA+scATAC) data, shared `REP_ID`, `PATH_to_barcodes` and `PATH_to_10x_ranger`.
2. multiple scRNA/scATAC/VISIUM replicates. One modality only.

```
SAMPLE	REP_ID	DATA_TYPE	PATH_to_barcodes	PATH_to_10x_bam	PATH_to_10x_ranger
HT001	U1	scRNA	barcodes.tsv.gz	s1.bam	s1/
HT001	U1	scATAC	barcodes.tsv.gz	s2.bam	s1/
```

### TODO
1. Detecting clonal LOH, from Numbat
```
In samples with high tumor purity (e.g., tumor cell lines) without matched normal cells, heterozygous SNPs are challenging to identify in regions of LoH, leading to decreased power of CNV detection. Regions of clonal LoH have decreased SNP density and can be identified by a specialized HMM (see detect_clonal_loh in function reference).
```
2. mapping resolutions, 2um to Xum?
3. remove genotype information in output file. Done


Rule of thumb from the docs: input/output/log/benchmark are relative to the working directory (--directory), while “other directives” (including script: and include:) are relative to the defining Snakefile.  

