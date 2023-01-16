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
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

print(individuals)

# construct paths for all possible outputs/targets, required for rule all
targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/06_ref_anno/" + s + "/sce_" + i + "-06"]
      targets = targets + [OUTPUT_BASE_PATH + "/reports/03_ref_annotation/" + s + "/ref_annotation_sample_report_" + i + ".html"]

targets = targets + [OUTPUT_BASE_PATH + "/reports/03_ref_annotation/ref_annotation_summary.html"]

#-------------------------------------------------------------------------------

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
  
"""
Preliminary cell type annotation using multiple reference datasets and
scmapcell and scmapcluster.
To get an overview of approx. cell type numbers and fraction contamination.
The longest sample takes around 1h for total 6 annotations.
"""
rule cell_type_annotation: 
    input: 
        sce_05 = OUTPUT_BASE_PATH + "/sce_objects/05_dimr/{species}/sce_{individual}-05",
        ref_baccin_sce = config["metadata"]["ref_baccin_sce"],
        ref_dahlin_sce = config["metadata"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["metadata"]["ref_dolgalev_sce"]
    output:
        sce_06 = OUTPUT_BASE_PATH + "/sce_objects/06_ref_anno/{species}/sce_{individual}-06"
    script:
        "scripts/06_annotation_scmap.R"  
      
#-------------------------------------------------------------------------------
# reports
  
rule make_sample_reports:
    input: 
        rules.cell_type_annotation.output
    output:
        OUTPUT_BASE_PATH + "/reports/03_ref_annotation/{species}/ref_annotation_sample_report_{individual}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "ref_annotation_sample_reports.Rmd" 

# generate all possible inputs for the summary
summary_inputs = []
for s in species:
  for i in individuals:
    if s in i:
      summary_inputs = summary_inputs+expand(rules.cell_type_annotation.output, species=s, individual=i)
      
rule make_summary:
    input:
        sce_06_pathlist = summary_inputs
    output:
        OUTPUT_BASE_PATH + "/reports/03_ref_annotation/ref_annotation_summary.html"
    params:
        samples_to_remove = config["samples_to_remove"],
        color_tables = TABLES_PATH,
        individuals = individuals,
        sce_functions = "../source/sce_functions.R" # we are in working dir
    script:
        #"scripts/testing_purposes.R"
        "ref_annotation_summary.Rmd" 
