#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import sys 

# paths from config
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["raw"])
print(OUTPUT_BASE_PATH)

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")
print(fractions) 
fractions = ["hsc"]

# construct paths for all possible outputs/targets, required for rule all
targets = []

for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_markers/markers_" + f]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_go/go_" + f]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_nmf/sce_nmf_" + f]
 
print(targets)

#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------

"""
get marker genes for each cluster and do go analysis on the marker genes.
Save both marker genes and GO results.
"""
rule go_analysis:
    input: 
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:
        results_markers = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_markers/markers_{fraction}",
        results_go = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_go/go_{fraction}"
    params:
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
    script:
        "scripts/11_go_analysis.R"
        
"""
Prepare NMF data and add to SCE objects for convenience.
"""
rule prepare_nmf:
    input: 
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_nmf/sce_nmf_{fraction}"
    script:
        "scripts/11_prepare_nmf.R"

