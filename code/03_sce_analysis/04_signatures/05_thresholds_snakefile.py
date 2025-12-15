#!/bin/python 

"""
TODO LATER
"""

#-------------------------------------------------------------------------------

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
METADATA_PATH = config["base"] + config["metadata_paths"]["metadata"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/05_thresholds"

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"]
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"]

#-------------------------------------------------------------------------------

targets = []
targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_p.rds"]
targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_t.rds"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

# recluster ts_hscs_progenitors using different sets of signatureso obtained
# after sliding the threshold on the nDGE analysis p value (by_p)
rule recl_ts_hscs_progenitors_by_p:
    input: 
        threshold_path = OUTPUT_DAT,
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_hscs_progenitors",
        ensembl_hum = ENSEMBL_HUM,
        ensembl_mus = ENSEMBL_MUS,
        metadata_path = METADATA_PATH
    output:
        seu_recl_output = OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_p.rds"
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        nr_cores = 10
    resources:
        mem_mb = 80000,
        queue = "medium-debian"
    threads: 10
    script:
        "05_thresholds/0X_recluster_by_p.R"

rule recl_ts_hscs_progenitors_by_t:
    input: 
        threshold_path = OUTPUT_DAT,
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_hscs_progenitors",
        ensembl_hum = ENSEMBL_HUM,
        ensembl_mus = ENSEMBL_MUS,
        metadata_path = METADATA_PATH
    output:
        seu_recl_output = OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_t.rds"
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        nr_cores = 10
    resources:
        mem_mb = 80000,
        queue = "medium-debian"
    threads: 10
    script:
        "05_thresholds/0X_recluster_by_t.R"