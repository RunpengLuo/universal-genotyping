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
