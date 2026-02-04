## Unified Genotyping Pipeline
```
mamba env create -f ./environment.yaml -p /path/to/envs/genotyping_env
conda activate /path/to/envs/genotyping_env

# run the pipeline
snakemake --cores <num_cores> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv

# generate snakemake performance report
snakemake --cores <num_cores> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv \
    --report report.html
```

### Overview
The universal genotyping pipeline supports genotyping, phasing, adaptive binning, RDR computation, etc., in either bulk (WGS/WES/lr-WGS) or single-cell (scRNA/ATAC-seq, Visium, Visium HD 3') sequencing data. it runs in two modes. For bulk mode, it produces necessary input files including BAF and RDR count matrices for running HATCHet. For single cell mode, it produces SNP/meta-SNP/bin-level allele count matrices for running CalicoST. Refer to `docs/bulk_mode.md` and `docs/single_cell_mode.md` for detailed configurations.

### Dependencies
Refer to `environment.yaml` for the complete list of required libraries. Additional phasing softwares:
1. Eagle2 `https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz`
2. Shapeit5 `https://github.com/odelaneau/shapeit5.git`, `https://github.com/odelaneau/shapeit5/releases`
3. LongPhase `https://github.com/twolinin/longphase`
