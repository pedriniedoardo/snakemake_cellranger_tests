# cellranger_test_01

* The samples are defined in the config file

## sample commands

```
snakemake --sdm conda -np all_test
snakemake --sdm conda -np all_cellranger
```
<br>

# cellranger_test_02

* The samples are defined in the config file
* better management of the output folder after cellranger
* added routine to move web summaries in one folder
* added routine to calculate the multiqc
* add pure intron-exon output

## sample commands

```
snakemake --sdm conda -np all_test
snakemake --sdm conda -np minimal_cellranger
snakemake --sdm conda -np default_cellranger
snakemake --sdm conda -np pure_IntronExon
```
<br>

# cellranger_test_03_multiruns

* The samples are defined in a `csv` file
* The pipeline can handle either single or multiple runs per sample.
* cellranger can be run in Single-run modality. Each sample/run combination generates a unique output
* cellranger can be run in Merged-run modality. Each sample generates a unique output. different runs are merged in one.

## Info
Because of the generation fo the `.mro` files during cellranger processing, it is necessary to add a unique ID to the sample name in the Single-run modality. At the moment I have decided to add a digit to the sample_name instead of the sample_run name. The link between digit and sample_run run can be seen in the web summary.
The multiQC run in Single-run modality is run over all the outputs once. The web summaries are divided by run.

## sample commands

### Single-run
```
snakemake --sdm conda -np cellranger_singlerun_default
```
### Merged-run
```
snakemake --sdm conda -np cellranger_multirun_default
snakemake --sdm conda -np multirun_pure_IntronExon
```