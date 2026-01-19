## Unified Single-cell/Visium Genotyping+Phasing+Pileup Workflow
```
snakemake -p --cores 1 --configfile config/config.yaml --directory <output>
```

### Dependencies
1. samtools
2. bcftools
3. cellsnp-lite
4. HetDetect `https://github.com/raphael-group/hetdetect.git`
5. Eagle2
6. Shapeit5

### Aux files
Required: reference panel SNPs VCF file and BAM files.
1. https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
    * used by Numbat, ~92MB
2. https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz
    * used by CalicoST, ~568MB

### Steps
If reference Het SNP file (genotyped from bulk) is provided, skip step 1-3:
1. pseudobulk (ignore barcode file) and genotype SNPs via cellsnp-lite mode1b.
2. select Het/Hom SNPs, via AF thresholding or HetDetect HMM
3. phase SNPs via population-based phasing methods
4. pile-up allele counts via cellsnp-lite mode1a.
5. postprocess allele count matrices.

### Sample file
Sample file has following modes:
1. multiome (GEX+ATAC) data, shared `PATH_to_barcodes` and `PATH_to_10x_ranger`.
2. multiple GEX/ATAC/VISIUM replicates. One modality only.

```
SAMPLE	REP_ID	DATA_TYPE	PATH_to_barcodes	PATH_to_10x_bam	PATH_to_10x_ranger
HT001	U1	GEX	barcodes.tsv.gz	s1.bam	s1/
HT001	U1	ATAC	barcodes.tsv.gz	s2.bam	s1/
```

### TODO
1. Detecting clonal LOH, from Numbat
```
In samples with high tumor purity (e.g., tumor cell lines) without matched normal cells, heterozygous SNPs are challenging to identify in regions of LoH, leading to decreased power of CNV detection. Regions of clonal LoH have decreased SNP density and can be identified by a specialized HMM (see detect_clonal_loh in function reference).
```
2. mapping resolutions, 2um to Xum?
3. remove genotype information in output file. Done


Rule of thumb from the docs: input/output/log/benchmark are relative to the working directory (--directory), while “other directives” (including script: and include:) are relative to the defining Snakefile.  

