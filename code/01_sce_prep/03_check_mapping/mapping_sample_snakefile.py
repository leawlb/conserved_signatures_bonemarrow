#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
OUTPUT_BASE = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["table"])
VALUES =  config["values"]

# specific report output path
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/01_sce_prep/"
OUTPUT_REP = OUTPUT_BASE + "/reports/01_sce_prep/03_check_mapping/"

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

# collect all paths for all possible outputs/targets, required for rule all
targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_DAT + "06_mapp/" + s + "/dmgs_" + i]
      targets = targets + [OUTPUT_REP + s + "/mapping_sample_report_" + i + ".html"]
      
targets = targets + [OUTPUT_DAT + "07_dmgs/dmgs_list_all"]

# local execution of non-demanding rules
localrules: all  

#-------------------------------------------------------------------------------

# define rules
rule all: # must contain all possible output paths from all rules below
    input:
        targets

#-------------------------------------------------------------------------------
# get the differentially mapped genes (dmgs) for each sample
rule get_sample_dmgs:
    input:
        sce_02 = OUTPUT_BASE + "/sce_objects/01_sce_prep/02_outl/{species}/sce_{individual}-02",
        sce_fg = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/fourgenomes/sce_objects/01_cellranger_output/{species}/sce_{individual}-01"
    output:
        dmgs = OUTPUT_DAT + "06_mapp/{species}/dmgs_{individual}"
    params:
        nr_hvgs = VALUES["preprocessing"]["nr_hvgs"],
        logFC_sample_mgs = VALUES["check_mapping"]["logFC_sample_mgs"],
        minPC_sample_mgs = VALUES["check_mapping"]["minPC_sample_mgs"]
    script:
        "scripts/06_sample_dmgs.R" 
    
# merge all dmgs into one list for exclusion during merging  
merge_dmgs_input = []
for s in species:
  for i in individuals:
    if s in i:
      merge_dmgs_input = merge_dmgs_input + expand(rules.get_sample_dmgs.output, species = s, individual = i)

rule merge_dmgs:
    input:
        dmgs = merge_dmgs_input
    output:
        dmg_list = OUTPUT_DAT + "07_dmgs/dmgs_list_all"
    script:
        "scripts/07_merge_dmgs.R"
       
rule make_sample_reports:
    input:
        sce_02 = OUTPUT_BASE + "/sce_objects/01_sce_prep/02_outl/{species}/sce_{individual}-02",
        sce_fg = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/fourgenomes/sce_objects/01_cellranger_output/{species}/sce_{individual}-01",
        dmgs = rules.get_sample_dmgs.output,
        dmg_list = rules.merge_dmgs.output
    output:
        OUTPUT_REP + "{species}/mapping_sample_report_{individual}.html"
    params:
        nr_hvgs = VALUES["preprocessing"]["nr_hvgs"],
        color_tables = TABLES_PATH
    script:
        "mapping_sample_report.Rmd" 
