#!/bin/python 

"""
Info in batch correction methods from:
Luecken et al. "Benchmarking atlas-level data integration in single-cell genomics", Nat Met 2022
Tran, Ang, Chevrier, Zhang et al. "A benchmark of batch effect correction methods for songle-cell RNA sequencing data", Genome Biologz 2020
"""
#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import sys 

# paths from config
OUTPUT_BASE_PATH = config["paths"]["output_dir"]

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
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

print(species)
# construct paths for all possible outputs/targets, required for rule all
targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_BASE_PATH + "/06_sglr/" + s + "/sce_" + i + "-06"]
      #targets = targets + [OUTPUT_BASE_PATH + "/07_mrge/sce_" + s + "-07"]
      #targets = targets + [OUTPUT_BASE_PATH + "/reports/03_annotation/annotation_species_report_" + s + ".html"]
  
if config["run_annotaion_summary"]:
  targets = targets + [OUTPUT_BASE_PATH + "/reports/03_annotation/annotation_summary.html"]
  
#-------------------------------------------------------------------------------

#localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
  
"""
Preliminary cell type annotation using multiple reference datasets and
SingleR.
To get an overview of approx. cell type numbers and fravtion contamination.
"""
rule cell_type_annotation: 
    input: 
        sce_04 = OUTPUT_BASE_PATH + "/04_norm/{species}/sce_{individual}-04"
    output:
        sce_06 = OUTPUT_BASE_PATH + "/06_sglr/{species}/sce_{individual}-06"
    params:
        ref_baccin_sce = config["metadata"]["ref_baccin_sce"],
        ref_dahlin_sce = config["metadata"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["metadata"]["ref_dolgalev_sce"],
        ref_lipka_sce = config["metadata"]["ref_lipka_sce"],
    script:
        "scripts/06_annotation_singleR.R"  
  
      
"""
Merge all datasets of one species into one species dataset
Input data is located in 02_preprocessing/04_norm/
"""
"""
rule merge_datasets_species:
    input: 
        sce_06 = OUTPUT_BASE_PATH + "/06_sglr/{species}/"
    output:
        sce_07 = OUTPUT_BASE_PATH + "/07_mrge/sce_{species}-07"
    params:
        individuals = individuals
    wildcard_constraints:
        species = "[a-z]+"
    script:
        "scripts/07_merge_datasets.R" 
"""  
