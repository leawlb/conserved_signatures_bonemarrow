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
individuals = get_list(metadata = METADATA, column = "Object_ID")

print(individuals)

# construct paths for all possible outputs/targets, required for rule all
targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/06_sglr/" + s + "/sce_" + i + "-06"]
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/06_sglr/" + s + "/preds/pred_all_" + i + "-06"]
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/06_sglr/" + s + "/preds/pred_hsc_" + i + "-06"]
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/06_sglr/" + s + "/preds/pred_str_" + i + "-06"]
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/03_reference_annotation/" + s + "/ref_annotation_sample_report_" + i + ".html"]

if config["run_ref_annotation_summary"]:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/03_reference_annotation/ref_annotation_summary.html"]

#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
  
"""
Preliminary cell type annotation using multiple reference datasets and
SingleR.
To get an overview of approx. cell type numbers and fraction contamination.
"""
rule cell_type_annotation: 
    input: 
        sce_05 = OUTPUT_BASE_PATH + "/sce_objects/05_dimr/{species}/sce_{individual}-05"
    output:
        sce_06 = OUTPUT_BASE_PATH + "/sce_objects/06_sglr/{species}/sce_{individual}-06",
        pred_hsc = OUTPUT_BASE_PATH + "/sce_objects/06_sglr/{species}/preds/pred_hsc_{individual}-06",
        pred_str = OUTPUT_BASE_PATH + "/sce_objects/06_sglr/{species}/preds/pred_str_{individual}-06",
        pred_all = OUTPUT_BASE_PATH + "/sce_objects/06_sglr/{species}/preds/pred_all_{individual}-06",
    params:
        ref_baccin_sce = config["metadata"]["ref_baccin_sce"],
        ref_dahlin_sce = config["metadata"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["metadata"]["ref_dolgalev_sce"]
    script:
        "scripts/06_annotation_singleR.R"  
      
#-------------------------------------------------------------------------------
# reports
     
rule make_ref_sample_reports:
    input: 
        rules.cell_type_annotation.output
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/03_reference_annotation/{species}/ref_annotation_sample_report_{individual}.html"
    params:
        nr_hvgs = config["metadata"]["values"]["nr_hvgs"],
        color_tables = TABLES_PATH
    script:
        "ref_annotation_sample_reports.Rmd" 
  
if config["run_ref_annotation_summary"]:
  print("run_ref_annotation_summary")      
rule make_ref_summary_report:
    input:
        sce_06_path = OUTPUT_BASE_PATH + "/sce_objects/06_sglr/",
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/03_reference_annotation/ref_annotation_summary.html"
    params:
        samples_to_remove = config["samples_to_remove"],
        color_tables = TABLES_PATH,
        individuals = individuals,
        species = species
    script:
        "ref_annotation_summary.Rmd" 
