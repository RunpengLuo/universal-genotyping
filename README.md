## Unified Genotyping Pipeline
This Snakemake pipeline takes a sample file listing data from same patient across different assays (bulkWGS (sr and lr), bulkWES, scRNA/ATAC-seq, VISIUM, VISIUM HD 3prime), a configuration file, and performs bulk genotyping, single-cell genotyping, or copy-typing preprocessing on request. It can do genotyping, phasing, adaptive binning, GC-corrected RDR computation, etc.,

### Quick-Start
```sh
mamba env create -f ./environment.yaml -p /path/to/envs/genotyping_env
conda activate /path/to/envs/genotyping_env

# dry-run to check input files are formatted properly
snakemake --cores 1 \
    --dry-run \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv

# run the pipeline
snakemake --cores <num_cores> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv

# generate report
snakemake --cores <num_cores> \
    -s /path/to/workflow/Snakefile \
    --configfile config/config.yaml \
    --directory <output> \
    --config sample_file=/path/to/sample_file.tsv \
    --report report.html
```

### Dependencies
Refer to `environment.yaml` for the complete list of dependecies libraries. Additional phasing softwares can be installed via mamba or GitHub epends on need:
1. Eagle2 `https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz`
2. Shapeit5 `https://github.com/odelaneau/shapeit5.git`, `https://github.com/odelaneau/shapeit5/releases`
3. LongPhase `https://github.com/twolinin/longphase`

### Configuration
To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file. An template sample file is provided as `config/samples.tsv`.

#### Parameters
TODO

