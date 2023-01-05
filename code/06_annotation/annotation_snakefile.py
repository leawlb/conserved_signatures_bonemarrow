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
  
species = get_list(metadata = METADATA, column = "Species_ID")
print(species) 

# construct paths for all possible outputs/targets, required for rule all
targets = [OUTPUT_BASE_PATH + "/sce_objects/14_scmp_cluster/"]

#for s in species:
  #targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/15_scmp_reference/sce_" + s + "-15"]

# summary reports
targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/06_annotation/annotation_summary_scmap.html"]

print(targets)

#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------


rule cluster_scmap_annotation:
    input: 
        sce_12_path = OUTPUT_BASE_PATH + "/sce_objects/12_lcls/"
    output:
        sce_14_path = OUTPUT_BASE_PATH + "/sce_objects/14_scmp_cluster/"
    params:
        species = species
    script:
        "scripts/14_scmap_mmus_cluster_annotation.R"
"""

rule reference_scmap_annotation:
    input: 
        sce_14 = rules.cluster_scmap_annotation.output + ["/sce_{species}-14"]
    output:
        sce_15 = OUTPUT_BASE_PATH + "/sce_objects/15_scmp_reference/sce_{species}-15"
    params:
        ref_baccin_sce = config["metadata"]["ref_baccin_sce"],
        ref_dahlin_sce = config["metadata"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["metadata"]["ref_dolgalev_sce"]
    script:
        "scripts/15_scmap_reference_celltype_annotation.R"
"""   

rule make_summary_report_scmap:
    input:
        sce_14_path = rules.cluster_scmap_annotation.output
    params:
       color_tables = TABLES_PATH,
       species = species
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/06_annotation/annotation_summary_scmap.html"
    script:
          "annotation_summary_scmap.Rmd" 
