#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/03_cci_analysis"
OUTPUT_REP = OUTPUT_BASE + "/cci_objects/reports/03_cci_analysis/04_cci_scores"

COLORS = config["base"] + config["metadata_paths"]["colors"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")

VALUES = config["values"]["05_cci_prep"]
LRDB = config["base"] + config["metadata_paths"]["lrdb_out"]

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_REP + "/scores_summary.html"]


#-------------------------------------------------------------------------------


localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

rule score_summary:
    input: 
        cci_met_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/10_ipis/ipi_{species}_{age}", species = species, age = age)
    output:
        OUTPUT_REP + "/scores_summary.html"
    params:
        species = species,
        age = age,
        colors_path = COLORS,
        colors = "../../source/colors.R"
    threads:
        4
    script:
        "cci_scores_summary.Rmd"  
