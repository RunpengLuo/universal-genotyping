## Unified WGS/WES/single-cell/Visium preprocessing

### Dependencies
1. samtools
2. bcftools
3. cellsnp-lite
4. HetDetect `https://github.com/raphael-group/hetdetect.git`

### Step 1. Genotype Het SNPs
Required: reference panel SNPs VCF file and BAM files.
1. https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
    * used by Numbat, ~92MB
2. https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz
    * used by CalicoST, ~568MB

Output: genotyped Het SNP VCF file

* For WGS/WES, genotype Het SNPs from matched-normal sample if available, otherwise genotype tumor sample via bcftools.
* For single-cell/Visium, pseudobulk (ignore barcode file) and genotype SNPs via cellsnp-lite.

* Infer Het SNPs via AF thresholding or HetDetect HMM
### Step 2. Reference/read-based phasing
Required: population panel BCF files, reference FASTA file.
* For long reads WGS, use long read phaser like LongPhase.
* For other modalities, use population-based phasing
    1. Eagle2
    2. SHAPEIT5
    3. Beagle5

    install numbat

### Observations
1. from Numbat
```
Detecting clonal LOH
In samples with high tumor purity (e.g., tumor cell lines) without matched normal cells, heterozygous SNPs are challenging to identify in regions of LoH, leading to decreased power of CNV detection. Regions of clonal LoH have decreased SNP density and can be identified by a specialized HMM (see detect_clonal_loh in function reference).
```



2um bam -> 16um
lose resolution at boundary?
# universal-genotyping
