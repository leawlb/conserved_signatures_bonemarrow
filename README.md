
This repository contains all scripts related to scRNAseq analysis of HSPCs
and niche cells from four distinct mouse species, starting from 
alignment with Cell Ranger up to figure creation.
Further, bulk RNAseq data from Neo1- and TrkB-stained stromal cells is 
pre-processed and analysed following experimental validation.


## 1. Configuration and Requirements

The general configuration is stored in `code/config.yaml`.
!!! Base paths MUST be adjusted before re-running code !!!

Minimum required base paths that need to be adjusted:
 - `code/config.yaml`
 - in `code/figures_rev`, check base_paths in all .Rmd files not run by the snakefile, and `OUTPUT_PATH` in the run_rmds.py file
 - in `metadata`, any .R scripts used strictly to deal with metadata (see below)
 - check `code/08_sce_brain/01_sce_brain/brain_snakefile.py`
 - check all files in `code/04_rare_celltypes` 

### Starting from raw scRNAseq files

For repeating alignment for scRNAseq data, required data are (more info below):

 - raw files (E-MTAB-15073) 
 - the N-masked reference genome (S-BSST2074) 
 - four species-specific reference genomes 

Additionally adjust base paths in:
 - `code/00_sce_cellranger/01_cellranger_main/config/config-interspecies-bonemarrow.yaml` 
 - `code/00_sce_cellranger/02_cellranger_fourgenomes/config/config-interspecies-bonemarrow.yaml` 
 - All paths in column "FastQ Path" in the `metadata/scRNAseq/00_alignment/metadata.csv` table must be adjusted to the appropriate paths after downloading the raw files.

For alignment of raw scRNAseq data with the N-masked reference genome or 
species-specific genomes navigate into
`code/00_sce_cellranger`, then into the appropriate subdirectories and 
see `README.txt` files there for set-up. 
These directories were adjusted from and added by Fritjof Lammers.

### Starting from processed scRNAseq files

The required data are (more info below):

 - processed files aligned with the N-masked reference genome (E-MTAB-15073) 
 
Navigate into `code/00_sce_cellranger/01_cellranger_main` in order to create 
SingleCellExperiment objects from the downloaded Cell Ranger output files.
See the `README.txt` file there for set-up.

Following steps in `code/01_sce_prep/01_preprocessing` include detection
and removal of genes differentially mapped between data aligned with the 
N-masked genome vs. data aligned with the species-specific genomes. 
See `code/01_sce_prep/01_preprocessing/preprocessing_snakefile.py` for more info.

These steps additionally require either (more info below):

 - four species-specific genomes to generate the required processed data via alignment
 - or the list of dmgs in `/metadata/scRNAseq/01_sce_prep/dmgs_list_all` that allows skipping the identification step. Copy this file into the appropriate folder designated in `code/01_sce_prep/01_preprocessing/preprocessing_snakefile.py` and adjust the snakefile to skip rules "get_sample_dmgs" and "make_dmg_reports".

Additionally adjust base paths in:
 - `code/00_sce_cellranger/01_cellranger_main/config/config-interspecies-bonemarrow.yaml`
 - in `code/00_sce_cellranger/01_cellranger_main/config/config-interspecies-bonemarrow.yaml` also adjust the cellranger_count path to appropriate path after download
 - `code/00_sce_cellranger/02_cellranger_fourgenomes/config/config-interspecies-bonemarrow.yaml` 

### Starting from annotated scRNAseq files

This is by far the easiest option for getting started. 

For starting from fully annotated objects (S-BSST2079) (more info below),
transfer the downloaded objects into the approriate folder analogous to 
`base` + `data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns` as 
determined in `code/config.yaml` and save them once using saveRDS but without 
".rds" file extension for downstream compatibility or adjust all affected
downstream paths as required.
This is because I mistakenly didn't use .rds extensions for a while.
Downstream analysis starts in folder 03.

### bulkRNAseq data

Start from raw bulkRNAseq fastq to repeat pre-processing including umi_tools,
cutadapt, alignment, featureCounts (`code/05_validation/01_bulkRNAseq`) etc.
Adjust `05_validation/01_bulkRNAseq/rules/01_link_files.smk` and the 
required paths in config as required.

Start from counts matrix and read into DESeq2 object for analysis via DESeq2
starting from `05_validation/01_bulkRNAseq/rules/10_deseq2.smk`.
Adjust data loading and DEseq2 object creation as required.


## 2. Snakemake set-up and execution 

### scRNAseq alignment

Navigate into `code/00_sce_cellranger`, then into the appropriate 
subdirectories and see `README.txt` files there for set-up. 

### scRNAseq analysis

For all steps starting from 01 install `snakemake_isbm.yml` micromamba:

```bash
micromamba env create -n snakemake_isbm -f snakemake_isbm.yml
```

Activate the environment: 

```bash
micromamba activate snakemake_isbm
```

Any other environments in `/env` will be automatically installed. 

Navigate to the appropriate folder (starting in 01_01) and run the snakemake
pipeline from that folder using the snakemake command specified 
in the `..._snakemake` txt files. Snakemake commands in the `..._snakemake` 
files are suited for the DKFZ cluster 
structure and must be adjusted to the local conditions before running.

Generally, follow the steps as indicated by numbers, even if some are missing 
(e.g. retired 05, 06, and 07 folders).
Run all analysis folders before running figure scripts.

### bulk RNAseq pre-processing and analysis

```bash
micromamba env create -n snakemake_bulk -f snakemake_bulk.yml
```

```bash
micromamba activate snakemake_bulk
```

Snakemake is set up differently in this directory, with a different version,
and relying on a snakemake profile, which has been adjusted to DKFZ cluster
requirements. 
Adjust snakemake setup as required.


## 2. Data

Data can be downloaded from BioStudies:

 - the N-masked reference genome, generated using this repository: https://github.com/fritjoflammers/snakemake-snpmasked-refgenome.git (https://doi.org/10.5281/zenodo.15516917) (S-BSST2074) 
 - raw data: scRNAseq fastq files (E-MTAB-15073)
 - processed data: Cell Ranger output files matrix.mtx, barcodes.tsv, features.tsv after alignment with the N-masked reference (E-MTAB-15073)
 - fully annotated data in .rds format containing cell type labels, normalised log-counts, and batch-corrected PC and UMAP coordinates (S-BSST2079)
 - raw data: bulkRNAseq fastq files (E-MTAB-17326)
 - processed data: bulkRNAseq counts matrix (featureCounts output) for all samples (E-MTAB-17326)
 

## 3. Metadata
 
Deposited public data used is also summarized in the Key resources table
of the associated publication.

Metadata required for running the code is either found in folder `metadata` 
(if they are relatively small, like assignment csvs) or can be generated there.
For big metadata files, .R scripts in the folder `metadata` that are required to obtain
this public data used for analysis.

 - `metadata/scRNAseq/00_alignment/metadata.csv`, which contains at least the same required information as the metadata table from E-MTAB-15073 but is formatted correctly for alignment and downstream analysis
 - cell type assignment lists
 - ensembl lists
 - gene lists
 - color schemes, etc.

The `metadata` folder should be copied in its entirety into the directory 
analogous to `base` + `data/` as determined in `code/config.yaml`. 


In addition, there is also some metadata that must be downloaded or generated manually.
 
- Four Cell Ranger reference genomes for species-specific alignment, generated from downloaded fasta and gtf files (http://ftp.ensembl.org/pub/release-94/) using Cell Ranger v3.1.0 mkref function. These can also be made available upon request, see associated publication:
  - GRCm38 (Ensembl release 94)
  - CAST_EiJ_v1 (Ensembl release 94)
  - SPRET_EiJ_v1 (Ensembl release 94)
  - CAROLI_EIJ_v1.1  (Ensembl release 94)
 
- Published datasets for reference annotation (see `/metadata/scRNAseq/01_sce_prep/references_raw/`)
  - Dahlin et al. (2018) 10X Genomics dataset was downloaded from the ABC portal (http://abc.sklehabc.com/)
  - Baccin et al. (2020): https://nicheview.shiny.embl.de
  - Dolgalev et al. (2021): https://osf.io/ne9vj/files/osfstorage 
 
- Published HSPC or Niche datasets for re-clustering, see `code/03_sce_analysis/04_signatures/00_prepare_datasets_snakefile.py` and `code/04_rare_celltypes` for more info. 

- The Bakken et al (2021) and Yao et al (2021) motor cortex datasets (`sample.combined_exc_4_species_integration.RDS`) must also be downloaded manually from https://data.nemoarchive.org/publication_release/Lein_2020_M1_study_analysis/Transcriptomics/sncell/10X/human/processed/analysis/analysis/M1/cross_species_integration/ 

- The mouse lemur dataset is available in figshare: https://figshare.com/projects/Tabula_Microcebus/112227

Third, even more other (meta)data may be downloaded automatically by running the code.

For example the entire github repository `mcclust` by Fritsch (2022) is 
downloaded as part of running the code so that its vi.dist function can be used.
See: https://github.com/cran/mcclust and 
`code/03_sce_analysis/04_signatures/01_reclustering_ow_snakefile.py`

Lastly, we thank Emmanouil Athanasiadis for sending his processed scRNAseq zebrafish HSPC data.

If you have any questions, please contact us (see corresponding information in the associated
publication).
