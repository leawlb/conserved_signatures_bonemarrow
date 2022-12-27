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
targets = [OUTPUT_BASE_PATH + "/sce_objects/10_hcls/sce_" + s + "-10" for s in species]
print(targets)

for s in species:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_scls/sce_" + s + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/05_manual_annotation/clustering/" + s + "/clustering_species_report1_" + s +".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/05_manual_annotation/clustering/" + s + "/clustering_species_report2_" + s +".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_scmp/sce_" + s + "-12"]

if config["run_manual_annotation_summary"]:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/05_manual_annotation/manual_annotation_species_summary.html"]
  
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------
# CLUSTERING (try two different methods)

rule hierarchical_clustering:
    input: 
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_{species}_Batch_exp_day-09"
    output:
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_hcls/sce_{species}-10"
    params:
        remove_mature_cells_clustering = config["remove_mature_cells_clustering"],
        number_k = config["values"]["clustering"]["number_k"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
    script:
        "scripts/10_hierarchical_clustering.R"
        
rule seurat_clustering:
    input: 
        sce_10 = rules.hierarchical_clustering.output
    output:
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_scls/sce_{species}-11"
    params:
        remove_mature_cells_clustering = config["remove_mature_cells_clustering"],
        separate_fractions_clustering = config["separate_fractions_clustering"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
    script:
        "scripts/11_seurat_clustering.R"
        
rule cluster_ref_annotation:
    input: 
        sce_11 = rules.seurat_clustering.output
    output:
        sce_12 = OUTPUT_BASE_PATH + "/sce_objects/12_scmp/sce_{species}-12"
    params:
        ref_baccin_sce = config["metadata"]["ref_baccin_sce"],
        ref_dahlin_sce = config["metadata"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["metadata"]["ref_dolgalev_sce"]
    script:
        "scripts/12_scmap_annotation.R"

rule make_clustering_species_report1:
    input:
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_{species}_Batch_exp_day-09",
        sce_10 = rules.hierarchical_clustering.output,
        sce_11 = rules.seurat_clustering.output,
        sce_12 = rules.cluster_ref_annotation.output
    params:
        color_tables = TABLES_PATH,
        number_k = config["values"]["clustering"]["number_k"]
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/05_manual_annotation/clustering/{species}/clustering_species_report1_{species}.html"
    script:
        "clustering_species_report_1.Rmd"
        
rule make_clustering_species_report2:
    input:
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_{species}_Batch_exp_day-09",
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/05_manual_annotation/clustering/{species}/clustering_species_report2_{species}.html"
    script:
        "clustering_species_report_2.Rmd"


#-------------------------------------------------------------------------------
# SUPPORT VECTOR MACHINES 


