#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import sys 

# paths and objects from config
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]
METADATA = pd.read_csv(config["metadata"]["raw"])

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

for f in fractions:
  targets = [OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_" + f + "-07"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/07_mrge_species/sce_" + s + "_" + f + "-07" for s in species]

#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
        
"""
Merge samples 
"""
rule merge_datasets_fractions:
    input: 
        sce_06 = expand(OUTPUT_BASE_PATH + "/sce_objects/06_ref_anno/{s}/", s = species)
    output:
        sce_07_hsc = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_hsc-07",
        sce_07_str = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_str-07"
    params:
        individuals = individuals,
        samples_to_remove = config["samples_to_remove"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
        species = species,
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    script:
        "scripts/07_merge_datasets_fractions.R" 
        
"""
Merge species 
"""
rule merge_datasets_species:
    input: 
        sce_06_path = OUTPUT_BASE_PATH + "/sce_objects/06_ref_anno/{species}/"
    output:
        sce_07_hsc_path = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_species/sce_{species}_hsc-07",
        sce_07_str_path = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_species/sce_{species}_str-07"
    params:
        individuals = individuals,
        samples_to_remove = config["samples_to_remove"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
        species = species,
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    script:
        "scripts/07_merge_datasets_species.R" 
        
        

