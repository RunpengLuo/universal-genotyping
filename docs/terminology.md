## cellsnp-lite definitions
1. in mode 1, cellsnp-lite takes BAM file and vcf file as input, and perform allele counting on REF and ALT allele over bi-allelic SNPs defined by the vcf file.
2. OTH is the counts that belongs to any allele besides REF and ALT allele defined by vcf file.
3. ALT allele count is AD, REF allele count is DP - AD. OTH+DP=total aligned read depth.

cases at true Het positions (based on sample)
1. sample has REF and ALT allele. OTH>0 indicates technical errors.
2. sample has REF and non-ALT allele. OTH>0 and AD=0, refers non-ALT allele counts, multi-allelic.
3. sample has non-REF and ALT allele. DP-AD=0, DP=AD, multi-allelic.
4. sample non-REF and non-ALT allele. DP=AD=0, filtered by cellsnp-lite.

cases at true Hom positions
1. 


I attempted to make a note for the discussion with Mike this afternoon and helped me to clarify some terms.

How cellsnp-lite and SNP genotyping works in current pipeline:
cellsnp-lite takes 1kgp vcf file (36.6M SNPs with minor allele frequency (MAF) > 0.0005), consider only bi-allelic REF/ALT SNPs defined by vcf file, and pileup reads aligned to REF (DP-AD) and ALT (named AD) alleles with barcodes ignored (pseudobulk). Any reads aligned to non-REF/ALT allele specified by 1kgp vcf file are counted as OTH . Any SNPs has pseudobulk DP<2 is excluded.
In Visium 3' ovarian dataset, there are 1202001 resulting SNPs, 54988 (~4%) of them have OTH >0.
OTH>0 cases exists, implies either the data has true non-REF/ALT allele, or due to technical errors like mis alignment or base calling error.
In both CalicoST and Numbat, a SNP is called Het if its minor allele-frequency MAF=min(AD, DP-AD)/DP >= 0.1 and minor allele count min(AD, DP-AD) >= 2. a SNP is called Hom if DP >=10 and min(AD, DP-AD)=0. Numbat also explicitly filters SNPs with OTH > 0, CalicoST didn’t take OTH into account, any nonzero OTH  SNPs can be silently assigned to Hom, Het, or rejected.
As the result, 28826=16279 (HET) + 12547 (HOM) SNPs are kept after genotyping step.
Mike suggests to validate Visium genotyped SNPs with paired WES genotyped SNPs, and potentially improve thresholding criteria or better genotyping strategy to achieve higher number of Het SNPs and also avoid to introduce many FP SNPs. (balance recall and precision).
From the literature besides thresholding strategy, SNP genotyping can be done by likelihood modeling on base calling error event and compute posteriors per genotype. There are also tools developed for single-cell SNV calling by making use of linkage disequilibrium from 1kgp haplotype blocks. Before trying new methods.
