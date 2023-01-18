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

# construct paths for all possible outputs/targets, required for rule all
targets = []
print(targets)

for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/10_hcls_fractions/sce_" + f + "-10"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/10_scls_fractions/sce_" + f + "-10"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_" + f + "-10"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/05_clustering_fractions/" + f + "/clustering_fraction_report_full_" + f +".html"]
  if config["calculate_k"]:
    targets = targets + [OUTPUT_BASE_PATH + "/reports/05_clustering_fractions/" + f + "/clustering_fraction_report_k_" + f +".html"]
    
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------

rule hierarchical_clustering:
    input: 
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnnc_fractions/sce_{fraction}_Object_ID-09"
    output:
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_hcls_fractions/sce_{fraction}-10"
    params:
        number_k = config["values"]["clustering"]["number_k"],
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        size_subs_hscs = config["values"]["clustering"]["size_subs_hscs"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
    script:
        "scripts/10_hierarchical_clustering.R"
        
rule seurat_clustering:
    input: 
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnnc_fractions/sce_{fraction}_Object_ID-09"
    output:
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/10_scls_fractions/sce_{fraction}-10"
    params:
        resolution = config["values"]["clustering"]["resolution"]
    script:
        "scripts/11_seurat_clustering.R"
        
rule louvain_clustering:
    input: 
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnnc_fractions/sce_{fraction}_Object_ID-09"
    output:
        sce_12 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    script:
        "scripts/12_louvain_clustering.R"

#-------------------------------------------------------------------------------   

if config["calculate_k"]:     
  rule make_report_k:
      input:
          sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnnc_fractions/sce_{fraction}_Object_ID-09",
      output:
          OUTPUT_BASE_PATH + "/reports/05_clustering_fractions/{fraction}/clustering_fraction_report_k_{fraction}.html"
      threads:
          20
      script:
          "clustering_fraction_report_k.Rmd"

rule make_report_full:
    input:
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnnc_fractions/sce_{fraction}_Object_ID-09",
        sce_10_h = rules.hierarchical_clustering.output,
        sce_10_s = rules.seurat_clustering.output,
        sce_10_l = rules.louvain_clustering.output
    params:
        color_tables = TABLES_PATH,
        size_subs_hscs = config["values"]["clustering"]["size_subs_hscs"],
        number_k = config["values"]["clustering"]["number_k"]
    output:
        OUTPUT_BASE_PATH + "/reports/05_clustering_fractions/{fraction}/clustering_fraction_report_full_{fraction}.html"
    script:
        "clustering_fraction_report_full.Rmd"
