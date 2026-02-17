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
Refer to `environment.yaml` for the complete list of dependecies libraries. Additional phasing softwares can be installed via mamba or GitHub depends on need:
1. Eagle2 `https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz`
2. Shapeit5 `https://github.com/odelaneau/shapeit5.git`, `https://github.com/odelaneau/shapeit5/releases`
3. LongPhase `https://github.com/twolinin/longphase`

### Usage
Refer to `docs/step-by-step.md` for step by step tutorial and `docs/manual.md` for detailed explanations of workflow rules. To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

### Sample file format
Sample file is a TSV-format file lists sample name, replicate IDs, and data paths. An template sample file is provided as `config/samples.tsv`. Columns and constraints are explained:
| Column  | Explanation |
| ------------- | ------------- |
| `SAMPLE`  | patient ID. required. All rows must have same value.  |
| `REP_ID` | replicate ID. required. `REP_ID` must be unique per row, excepted for multiome (scRNA+scATAC) dataset.  |
| `assay_type` | Assay type. required. Choices: {`bulkWGS`, `bulkWES`, `scATAC`, `scRNA`, `VISIUM`, `VISIUM3prime`} |
| `sample_type` | Sample type. required. Choices: {`normal`, `tumor`} |
| `PATH_to_bam` | Path to .bam file. required. |
| `PATH_to_barcodes` | Path to 10x `barcodes.tsv.gz` file. required for non-bulk data. Can be found under `outs/filtered_feature_bc_matrix/`.|
| `PATH_to_10x_ranger` | Path to 10x cell-ranger or space-ranger `outs` directory. required for non-bulk data. `filtered_feature_bc_matrix*` and `{atac_}*fragments.tsv.gz` (scATAC) must present under `<PATH_to_10x_ranger>/`. |
| `PATH_to_ref_annotations` | Path to reference annotation TSV-format file. columns `BARCODE` and `<ref_label>` are required. |


