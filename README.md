
This repository contains all scripts related to scRNAseq analysis starting from 
alignment with Cell Ranger up to figure creation for this publication: 

DOI


### 1. Configuration

The configuration is stored in config.yaml.
Base paths MUST be adjusted before re-running code. 
Adjust paths in: 

 - `config-interspecies-bonemarrow.yaml` (2x)

 - `config.yaml`

 - all .Rmd files in the `figures` folder 
 

### 2. Alignment

For alignment of raw data with the N-masked reference genome (00_01) or 
species-specific genomes (00_02) navigate into the appropriate directories 
and see README.txt files there. 
These directories were adjusted from and added by Fritjof Lammers.


### 3. Snakemake set-up and execution

For all steps starting from 01 install snakemake_isbm.yml micromamba:

`micromamba env create -n snakemake_isbm -f snakemake_isbm.yml`

Activate the environment: 

`micromamba activate snakemake_isbm`

Navigate to the appropriate folder (starting in 01_01) and run the snakemake
pipeline from that folder using the snakemake command specified 
in the `..._snakemake` txt files. 

If no similar cluster LSF structure to the DKFZ cluster is available,
snakemake commands in the `..._snakemake` files must be adjusted to the 
local conditions. 


### 4. Data

Data can be downloaded from ArrayExpression or BioStudies:

 - Raw data: fastq files (E-MTAB-15073)
 
 - Processed data: matrix.mtx, barcodes.tsv, features.tsv after alignment (E-MTAB-15073)
 
 - N-masked reference genome, generated using this repository: LINK () #TODO: INSERT LINKS 
 
 - Fully annotated data in .Rds format containing cell type labels, normalised log-counts, and batch-corrected PC and UMAP coordinates (LINK) #TODO: INSERT LINK 
 
 
### 5. Metadata
 
Metadata used for running the code is in folder `metadata`:

 - `metadata.csv` contains at least the same required information as the metadata table from E-MTAB-15073 but is formatted correctly for alignment and downstream analysis (see “/metadata/scRNAseq/00_alignment”)
 
 - cell type assignment lists
 
 - gene lists
 
 - color schemes, etc.

The `metadata` folder should be copied in its entirety into the directory used to store data as determined in `config.yaml` for smooth running. 
 

Other (meta)data must be downloaded manually or can be made available upon 
request to l.woelbert[at]dkfz-heidelberg.de:
 
Four Cell Ranger reference genomes for species-specific alignment, 
generated from downloaded fasta and gtf files 
(http://ftp.ensembl.org/pub/release-94/) using Cell Ranger v3.1.0 mkref function
 - GRCm38 (Ensembl release 94)
 - CAST_EiJ_v1 (Ensembl release 94)
 - SPRET_EiJ_v1 (Ensembl release 94)
 - CAROLI_EIJ_v1.1  (Ensembl release 94)
 
Published datasets for reference annotation 
(see “/metadata/scRNAseq/01_sce_prep/references_raw/”)
 - Dahlin et al. (2018) 10X Genomics dataset was downloaded from the ABC portal (http://abc.sklehabc.com/)
 - Baccin et al. (2020): https://nicheview.shiny.embl.de
 - Dolgalev et al. (2021): https://osf.io/ne9vj/files/osfstorage 
 
Published datasets for re-clustering, instructions in 03/04/00 Snakefile
 
Other metadata may be downloaded automatically by running the code.
