## Unified genotyping pipeline strategies
This unified pipeline produces allelic and total data matrices for bulk and single-cell CNA callers, including HATCHet3, Copy-typing, CalicoST, and cnaster.

### Support HATCHet3
Input: 1 matched-normal sample and multiple tumor sample, short-read WGS/WES, or long-read WGS data.
1. `genotype_snps_bulk` takes matched-normal BAM file and SNP locations from SNP panel/database, outputs per-chromosome bi-allelic genotyped SNPs via bcftools.
2. `phase_snps_<phaser>` takes bi-allelic genotyped SNPs and phasing panel (or matched-normal BAM file for long reads), output phased Het SNPs after merging&filtering.
3. `pileup_snps_cellsnp_lite_bulk` takes phased Het SNPs and BAM files, output per-sample SNP-level DP (ALT+REF) and AD(ALT) allele count vectors via cellsnp-lite.
4. `postprocess_bulk` takes phased Het SNPs and allele count vectors, and additional BED files including following:
    1. `region_bed`: BED file(s) contains non-overlapping valid intervals (e.g., non-centromeric regions).
    2. `block_bed`: BED file contains non-overlapping aggregation unit blocks like gene blocks, aggregated bins cannot contain partial blocks.

    And do the following:
    1. concat per-sample count vectors into SNP by sample count matrices.
    2. filter SNPs based on `region_bed`, per-sample min-depth, beta-posterior credible interval AF test on normal sample. If `block_bed` is given, filter SNPs not in `block_bed`.
    3. convert to A/B phased matrices.
    4. set tag `region_id` and assign per-SNP sub-intervals based on adjacent midpoints along chromosome.
    5. TODO handle `block_bed`.
    6. assign `binom_id` based on pairwise AF binom test. same id allows to aggregate.
    6. `PS` given by phaser or estimated by switch probs inferred from genetic map with threshold `max_swithprob`.
    7. group SNPs by (`region_id`, `binom_id`, `PS`), per group, assign `meta_id` and `bbc_id` refers to n=2 SNP blocks and SNP blocks with MSR and MSPB criteria.
    8. build segment df `segments0` and `segments`, record left and right most SNP positions as `POS0L` and `POS0R`, expand to segment BED-like intervals `START`, `END` either based on adjacent midpoint (expanded regions, HATCHet convention) or set `START=POS0L` and `END=POS0R+1`.
    9. aggregate counts into segment by sample phased matrices according to step 4.8.
5. `run_mosdepth_bulk` takes `segments` BED file and computes per-segment read depth using BAM files.
6. `count_reads_bulk` convert read depth BED files into segment by sample read depth matrix.
7. `compute_rdr_bulk` computes read-depth ratio, and perform gc normalization (optional).

Output:
    1. 
