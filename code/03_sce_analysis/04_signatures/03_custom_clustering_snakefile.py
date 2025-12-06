#!/bin/python 

"""
This script is to perform custom clustering for three of the eight published
datasets used to evaluate the signature genes.

The three chosen datasets are ts_bone_marrow, mus_weinreb_hspc, and 
ts_all_stromal.

They were chosen because they have relatively low resolution of cell types/
clusters after I removed cells not found in our datase.
For example, ts_bone_marrow contained all bone marrow cells, including mature
cells, which I had however removed, so the remaining progenitors, originally
making up only a small fraction, are labelled in low resolution.

To fairly compare to our reclustering with optimized clustering resolution, 
I will perform a custom clustering and cell type annotation based on my 
expertise in computational analysis of hematopoietic and niche cells.

Because it is much simpler and each dataset needs to be treated highly
individually, I will perform the analysis in .Rmd scripts in contrast to
the usual .R scripts, one for each dataset.
"""

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/03_custom"
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/03_custom_clustering/"

COLORS = config["base"] + config["metadata_paths"]["colors"]

METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_REP + "/custom_ts_bonemarrow.html"]


#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# tabula sapiens adult bone marrow

rule custom_ts_bonemarrow:
    input: 
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_bone_marrow",
        ensembl_hum = config["base"] + config["metadata_paths"]["ensembl_hum"]
    output:
        OUTPUT_REP + "/custom_ts_bonemarrow.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "03_scripts_custom/custom_ts_bonemarrow.Rmd"
