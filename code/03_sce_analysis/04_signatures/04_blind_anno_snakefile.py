#!/bin/python 

"""
This script is to perform blind annotation for the selected other datasets
previously clustered using only signature genes.

The chosen datasets are ts_bone_marrow, ts_hscs_progenitors, mus_tik_stromal
and li_all_stromal.

Because it is much simpler and each dataset needs to be treated highly
individually, I will perform the analysis in .Rmd scripts in contrast to
the usual .R scripts, one for each dataset.
Then, I will save the object from within the .Rmd script.
I am also loading the required datasets within each separate .Rmd script.
"""

#-------------------------------------------------------------------------------

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
METADATA_PATH = config["base"] + config["metadata_paths"]["metadata"]
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/04_blind_anno"

"""
METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
"""

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_REP + "/blind_anno_ts_bonemarrow.html"]
targets = targets + [OUTPUT_REP + "/blind_anno_ts_hscs_progenitors.html"]
targets = targets + [OUTPUT_REP + "/blind_anno_mus_tik_stromal.html"]
targets = targets + [OUTPUT_REP + "/blind_anno_li_all_stromal.html"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

# tabula sapiens adult bone marrow
rule blind_ts_bonemarrow:
    input: 
        base_path = OUTPUT_BASE,
        metadata_path = METADATA_PATH
    output:
        OUTPUT_REP + "/blind_anno_ts_bonemarrow.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    script:
        "04_blind_annotation/blind_anno_ts_bonemarrow.Rmd"

# tabula sapiens fetal hspc 
rule blind_ts_hscs_progenitors:
    input: 
        base_path = OUTPUT_BASE,
        metadata_path = METADATA_PATH
    output:
        OUTPUT_REP + "/blind_anno_ts_hscs_progenitors.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    script:
        "04_blind_annotation/blind_anno_ts_hscs_progenitors.Rmd"

# Tikhonova et al mouse stromal cells
rule blind_mus_tik_stromal:
    input: 
        base_path = OUTPUT_BASE,
        metadata_path = METADATA_PATH
    output:
        OUTPUT_REP + "/blind_anno_mus_tik_stromal.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    script:
        "04_blind_annotation/blind_anno_mus_tik_stromal.Rmd"

# Li et al stromal cells
rule blind_li_all_stromal:
    input: 
        base_path = OUTPUT_BASE,
        metadata_path = METADATA_PATH
    output:
        OUTPUT_REP + "/blind_anno_li_all_stromal.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    script:
        "04_blind_annotation/blind_anno_li_all_stromal.Rmd"
