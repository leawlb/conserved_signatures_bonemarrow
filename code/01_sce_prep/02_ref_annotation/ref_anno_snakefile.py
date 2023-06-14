#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
OUTPUT_BASE = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = "/sce_objects/01_sce_prep/"
OUTPUT_REP = "/reports/01_sce_prep/02_ref_anno/"

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

#-------------------------------------------------------------------------------

# construct paths for all possible outputs/targets, required for rule all
targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_BASE + OUTPUT_DAT + "05_rfan/" + s + "/sce_" + i + "-05"]
      targets = targets + [OUTPUT_BASE + OUTPUT_REP +  s + "/02_ref_anno_sample_report_" + i + ".html"]

targets = targets + [OUTPUT_BASE + OUTPUT_REP + "ref_anno_summary.html"]

#-------------------------------------------------------------------------------

localrules: all  

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
rule sample_anno: 
    input: 
        sce_input = OUTPUT_BASE + OUTPUT_DAT + "04_dimr/{species}/sce_{individual}-04",
        ref_baccin_sce = config["metadata"]["ref_baccin_sce"],
        ref_dahlin_sce = config["metadata"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["metadata"]["ref_dolgalev_sce"]
    output:
        sce_05 = OUTPUT_BASE + OUTPUT_DAT + "05_rfan/{species}/sce_{individual}-05"
    script:
        "scripts/05_ref_anno_scmap.R"  
      
#-------------------------------------------------------------------------------
# reports
  
rule make_sample_reports:
    input: 
        rules.sample_anno.output
    output:
        OUTPUT_BASE + OUTPUT_REP + "{species}/02_ref_anno_sample_report_{individual}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "ref_anno_sample_reports.Rmd" 

summary_inputs = []
for s in species:
  for i in individuals:
    if s in i:
      summary_inputs = summary_inputs+expand(rules.sample_anno.output, species=s, individual=i)
      
rule make_summary:
    input:
        sce_05_pathlist = summary_inputs
    output:
        OUTPUT_BASE + OUTPUT_REP + "ref_anno_summary.html"
    params:
        samples_to_remove = config["samples_to_remove"],
        color_tables = TABLES_PATH,
        individuals = individuals
    script:
        "ref_anno_summary.Rmd" 
