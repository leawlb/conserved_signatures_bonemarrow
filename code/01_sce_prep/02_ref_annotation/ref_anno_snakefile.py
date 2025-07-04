#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/01_sce_prep"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/01_sce_prep/02_ref_anno"

print(OUTPUT_BASE)

METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
print(COLORS_REF)

#-------------------------------------------------------------------------------

targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_DAT + "/07_rfan/" + s + "/sce_" + i + "-07"]
      targets = targets + [OUTPUT_REP + "/" + s + "/02_ref_anno_report_" + i + ".html"]

targets = targets + [OUTPUT_REP + "/ref_anno_summary.html"]
#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
  
"""
Preliminary cell type annotation using multiple reference datasets and
scmapcell and scmapcluster.
To get an overview of approx. cell type numbers.
The longest sample takes around 1h for total 6 annotations.
"""
rule ref_anno: 
    input: 
        sce_input = OUTPUT_DAT + "/06_dimr/{species}/sce_{individual}-06",
        ref_baccin_sce = config["base"] + config["metadata_paths"]["ref_baccin_sce"],
        ref_dahlin_sce = config["base"] + config["metadata_paths"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["base"] + config["metadata_paths"]["ref_dolgalev_sce"]
    resources:
        mem_mb=50000
    output:
        sce_output = OUTPUT_DAT + "/07_rfan/{species}/sce_{individual}-07"
    script:
        "scripts/07_ref_anno_scmap.R"  
      
#-------------------------------------------------------------------------------
rule make_sample_reports:
    input: 
        rules.ref_anno.output
    resources:
        mem_mb=5000
    output:
        OUTPUT_REP + "/{species}/02_ref_anno_report_{individual}.html"
    params:
        colors_ref_path = COLORS_REF,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R",
        functions = "../../source/sce_functions.R"
    script:
        "ref_anno_sample_reports.Rmd" 

summary_inputs = []
for s in species:
  for i in individuals:
    if s in i:
      summary_inputs = summary_inputs+expand(rules.ref_anno.output, species=s, individual=i)
      
rule make_summary:
    input:
        sce_input_pathlist = summary_inputs
    resources:
        mem_mb=20000
    output:
        OUTPUT_REP + "/ref_anno_summary.html"
    params:
        samples_to_remove = config["samples_to_remove"],
        plotting = "../../source/plotting.R",
        individuals = individuals
    script:
        "ref_anno_summary.Rmd" 
