# snakemake-cellranger


This is the Cellranger workflow for the Odomlab. 

![rule graph](docs/rulegraph.png)


The pipeline and this documentation is work in progress and incomplete. Feel free to improve it by 
[forking the repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks) and making 
a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork). 

## Features

 - Read OTP-exported metadata to gather sample information (OTP/DKFZ format).
 - Link FASTQ files as CellRanger compatible filenames
 - Run CellRanger count for Gene Expression datasets.   
 - Store (subsetted) metadata CSV for downstream use. 
 - generate SingleCellExperiment objects per sample

## Feature ideas

- Add Cellranger ATAC
- Add Spaceranger
- Create a Seurat object and store as Rdata file.
- Create a SingleCellExperiment object and store as Rdata file.
- Run MultiQC on available FastQC file.
- *...*


## How to get

Clone this repository.

On the DKFZ cluster, you need to use HTTPS because SSH connections are blocked by the proxy server. 

To avoid typing your Github password each time, you may use Github's [Personal Access tokens](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) and store it in your 
Shell. 

## How to configure

An example configuration is in `config/config-interspecies_bonemarrow.yaml`.

To adapt the workflow, make a copy of this file and adapt it to your needs.

Minimal changes needed are: 

 * Metadata: Path to a OTP-exported metadata in `["metadata"]["table"]`. Multiple files can be listed. 
 * Output directory: `["paths"]["output_dir"]` - where the results should be stored
 * references: Choose directory where the N-merged reference genome is stored
 * cellranger_path path: install Cellranger 6.1.1 and store path

## How to run 

Use snakemake_cellranger.yaml to create an environment with required packages
for running snakemake.

```bash
micromamba create -f snakemake_cellranger.yaml
```

Set channel priority to strict.

```bash
conda config --set channel_priority strict
```

You may call the pipeline as follows in the directory where you cloned it. 
Adjusted to the DKFZ cluster architecture, may need to adjust as required.

```bash
snakemake --cluster "bsub -n16 -q long-debian -R rusage[mem=200GB]" -p -j4 -c42 --configfile config/config-interspecies-bonwmarrow.yaml --use-conda --conda-frontend conda
```

 - `--cluster` may change depending on the computational footprint of your analyses
 - `--configfile` should point to your personal configuration
 
Or directly adjusted to the data and DKFZ cluster architecture, may need to alter `workflow/Snakefile` as required:

```bash
snakemake --cluster "bsub -n{threads} -q {resources.queue} -R rusage[mem={resources.mem_mb}]" -p -j66 -c42 --configfile ../config/config-interspecies-bonemarrow.yaml --latency-wait 60 --use-conda --conda-frontend conda -n
```
