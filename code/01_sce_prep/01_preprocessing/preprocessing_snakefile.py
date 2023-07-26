#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths and objects from config
OUTPUT_BASE = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

METADATA = pd.read_csv(config["metadata"]["table"])
VALUES = config["values"]["preprocessing"]

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/01_sce_prep/"
OUTPUT_REP = OUTPUT_BASE + "/reports/01_sce_prep/01_preprocessing/"

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

#-------------------------------------------------------------------------------

# collect all paths for required outputs/targets, required for rule all
targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_DAT + "01_drop/" + s + "/sce_" + i + "-01"]
      targets = targets + [OUTPUT_DAT + "02_mapp/" + s + "/dmgs_" + i]
      targets = targets + [OUTPUT_DAT + "04_outl/" + s + "/sce_" + i + "-04"]
      targets = targets + [OUTPUT_DAT + "05_norm/" + s + "/sce_" + i + "-05"]
      targets = targets + [OUTPUT_DAT + "06_dimr/" + s + "/sce_" + i + "-06"]
      targets = targets + [OUTPUT_REP + "qc/" + s + "/preprocessing_qc_report_" + i + ".html"]
     # targets = targets + [OUTPUT_REP + "dmgs/" + s + "/preprocessing_dmg_report_" + i + ".html"]

#targets = targets + [OUTPUT_DAT + "03_dmgs/dmgs_list_all"]
      
if config["run_preprocessing_summary"]:
  targets = targets + [OUTPUT_REP + "qc/preprocessing_qc_summary.html"]

# local execution of non-demanding rules
localrules: all  

#-------------------------------------------------------------------------------
# define rules
rule all: # must contain all possible output paths from all rules below
    input:
        targets
        
#-------------------------------------------------------------------------------
# remove empty droplets and doublets to retain only single cells
rule remove_droplets:
    input: 
        sce_input = OUTPUT_BASE + config["paths"]["cellranger_output"] + "/{species}/sce_{individual}-01"
    output:
        sce_output = OUTPUT_DAT + "01_drop/{species}/sce_{individual}-01"
    params:
        cutoff_umis = VALUES["cutoff_umis"],
        cutoff_doublets = VALUES["cutoff_doublets"]
    script:
        "scripts/01_remove_droplets.R"

# get the differentially mapped genes (dmgs) for each sample 
rule get_sample_dmgs:
    input:
        sce_input = rules.remove_droplets.output,
        sce_fg = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/fourgenomes/sce_objects/01_cellranger_output/{species}/sce_{individual}-01"
    output:
        dmgs = OUTPUT_DAT + "02_mapp/{species}/dmgs_{individual}"
    params:
        nr_hvgs = VALUES["nr_hvgs"],
        logFC_sample_dmgs = VALUES["logFC_sample_dmgs"],
        minPC_sample_dmgs = VALUES["minPC_sample_dmgs"]
    script:
        "scripts/02_sample_dmgs.R" 
    
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
        dmg_list = OUTPUT_DAT + "03_dmgs/dmgs_list_all"
    script:
        "scripts/03_merge_dmgs.R"
        
# find outlier cells with low quality and remove them, also remove dmgs 
rule remove_outliers_dmgs:
    input: 
        sce_input = rules.remove_droplets.output,
        dmg_list = rules.merge_dmgs.output
    output:
        sce_output = OUTPUT_DAT + "04_outl/{species}/sce_{individual}-04"
    params:
        cutoff_sum = VALUES["cutoff_sum"],
        cutoff_detected = VALUES["cutoff_detected"],
        cutoff_mitos = VALUES["cutoff_mitos"]
    script:
        "scripts/04_remove_outliers_dmgs.R"
        
# normalize expression on sample level (for NMFs and Annotation)
rule normalize_expr:
    input: 
        sce_input = rules.remove_outliers_dmgs.output
    output:
        sce_output = OUTPUT_DAT + "05_norm/{species}/sce_{individual}-05"
    script:
        "scripts/05_normalize_expr.R"
        
# extract hvgs and reduce dimensions for QC visualisation on sample level 
rule reduce_dims:
    input: 
        sce_input = rules.normalize_expr.output
    output:
        sce_output = OUTPUT_DAT + "06_dimr/{species}/sce_{individual}-06" 
    params:
        nr_hvgs = VALUES["nr_hvgs"]
    script:
        "scripts/06_reduce_dims.R"
        
#-------------------------------------------------------------------------------
# REPORTS

# report file for dmgs 
rule make_dmg_reports:
    input:
        sce_input = rules.remove_droplets.output,
        sce_fg = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/fourgenomes/sce_objects/01_cellranger_output/{species}/sce_{individual}-01",
        dmgs = rules.get_sample_dmgs.output,
        dmg_list = rules.merge_dmgs.output
    output:
        OUTPUT_REP + "dmgs/{species}/preprocessing_dmg_report_{individual}.html"
    params:
        nr_hvgs = VALUES["nr_hvgs"],
        color_tables = TABLES_PATH
    script:
        "preprocessing_dmg_reports.Rmd" 

# report file for QC
rule make_qc_reports:
    input:
        sce_input = OUTPUT_BASE + config["paths"]["cellranger_output"] + "/{species}/sce_{individual}-01",
        sce_drop = rules.remove_droplets.output,
        sce_outl = rules.remove_outliers_dmgs.output,
        sce_norm = rules.normalize_expr.output,
        sce_dimr = rules.reduce_dims.output
    output:
        OUTPUT_REP + "qc/{species}/preprocessing_qc_report_{individual}.html"
    params:
        cutoff_umis = VALUES["cutoff_umis"],
        cutoff_doublets = VALUES["cutoff_doublets"],
        cutoff_sum = VALUES["cutoff_sum"],
        cutoff_detected = VALUES["cutoff_detected"],
        cutoff_mitos = VALUES["cutoff_mitos"],
        nr_hvgs = VALUES["nr_hvgs"],
        color_tables = TABLES_PATH
    script:
        # .Rmd files should be stored in snakefile working directory
        "preprocessing_qc_reports.Rmd" 
        
"""
Make one summary report on all files
Because this rule is dependend on previous outputs but this is not reflected
in input, it's generally turned off in config to avoid premature job submission
"""
print(config["run_preprocessing_summary"])
if config["run_preprocessing_summary"]:
  rule make_summary_report:
      input:
          sce_input_path = OUTPUT_DAT + "01_drop/"
      output:
          OUTPUT_REP + "qc/preprocessing_qc_summary.html"
      params:
          cutoff_sum = VALUES["cutoff_sum"],
          cutoff_detected = VALUES["cutoff_detected"],
          cutoff_mitos = VALUES["cutoff_mitos"],
          individuals = individuals,
          metadata = config["metadata"]["table"],
          color_tables = TABLES_PATH
      script:
          "preprocessing_qc_summary.Rmd" 
