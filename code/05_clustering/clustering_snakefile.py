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
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_lcls/sce_" + s + "-12"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_svms/objects_hierarchical/objects_" + s]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_svms/objects_seurat/objects_" + s]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_svms/objects_louvain/objects_" + s]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/" + s + "/clustering_species_report_k_" + s +".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/" + s + "/clustering_species_report_full_" + s +".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/" + s + "/clustering_species_report_svm_" + s +".html"]

if config["run_manual_annotation_summary"]:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/clustering_summary.html"]
  
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------

rule hierarchical_clustering:
    input: 
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_{species}_Batch_exp_day-09"
    output:
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_hcls/sce_{species}-10"
    params:
        remove_mature_cells_clustering = config["remove_mature_cells_clustering"],
        number_k = config["values"]["clustering"]["number_k"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
        separate_fractions_clustering = config["separate_fractions_clustering"]
    script:
        "scripts/10_hierarchical_clustering.R"
        
rule seurat_clustering:
    input: 
        sce_10 = rules.hierarchical_clustering.output
    output:
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_scls/sce_{species}-11"
    params:
        separate_fractions_clustering = config["separate_fractions_clustering"]
    script:
        "scripts/11_seurat_clustering.R"
        
rule louvain_clustering:
    input: 
        sce_11 = rules.seurat_clustering.output
    output:
        sce_12 = OUTPUT_BASE_PATH + "/sce_objects/12_lcls/sce_{species}-12"
    params:
        separate_fractions_clustering = config["separate_fractions_clustering"]
    script:
        "scripts/12_louvain_clustering.R"
   
rule run_svm_all:
    input: 
        sce_12 = rules.louvain_clustering.output
    output:
        sce_test = OUTPUT_BASE_PATH + "/sce_objects/13_svms/sce_test/sce_test_{species}",
        objects_hierarchical = OUTPUT_BASE_PATH + "/sce_objects/13_svms/objects_hierarchical/objects_{species}",
        objects_seurat = OUTPUT_BASE_PATH + "/sce_objects/13_svms/objects_seurat/objects_{species}",
        objects_louvain = OUTPUT_BASE_PATH + "/sce_objects/13_svms/objects_louvain/objects_{species}"
    params:
        use_kernel = config["values"]["clustering"]["use_svm_kernel"]
    threads:
        20
    script:
        "scripts/13_svm.R"
        
#-------------------------------------------------------------------------------   
        
rule make_clustering_species_report_k:
    input:
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_{species}_Batch_exp_day-09",
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/{species}/clustering_species_report_k_{species}.html"
    threads:
        20
    script:
        "clustering_species_report_k.Rmd"

rule make_clustering_species_report_full:
    input:
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_{species}_Batch_exp_day-09",
        sce_12 = rules.louvain_clustering.output
    params:
        color_tables = TABLES_PATH,
        number_k = config["values"]["clustering"]["number_k"]
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/{species}/clustering_species_report_full_{species}.html"
    script:
        "clustering_species_report_full.Rmd"
  
rule make_clustering_species_report_svms:
    input:
        sce_test = rules.run_svm_all.output.sce_test,
        objects_hierarchical = rules.run_svm_all.output.objects_hierarchical,
        objects_seurat = rules.run_svm_all.output.objects_seurat,
        objects_louvain = rules.run_svm_all.output.objects_louvain
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/{species}/clustering_species_report_svm_{species}.html"
    script:
        "clustering_species_report_svm.Rmd"      

if config["run_clustering_summary"]:
  rule make_summary_report:
      input:
          sce_12_path = OUTPUT_BASE_PATH + "/sce_objects/12_lcls",
          objects_13_path = OUTPUT_BASE_PATH + "/sce_objects/13_svms"
      params:
          species = species
      output:
          OUTPUT_BASE_PATH + "/sce_objects/reports/05_clustering/clustering_summary.html"
      script:
          "clustering_summary.Rmd" 
