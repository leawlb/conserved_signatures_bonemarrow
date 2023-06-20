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
      targets = targets + [OUTPUT_DAT + "02_outl/" + s + "/sce_" + i + "-02"]
      targets = targets + [OUTPUT_DAT + "03_norm/" + s + "/sce_" + i + "-03"]
      targets = targets + [OUTPUT_DAT + "04_dimr/" + s + "/sce_" + i + "-04"]
      targets = targets + [OUTPUT_REP + s + "/preprocessing_sample_report_" + i + ".html"]
      
if config["run_preprocessing_summary"]:
  targets = targets + [OUTPUT_REP + "preprocessing_summary.html"]

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
        sce_01 = OUTPUT_DAT + "01_drop/{species}/sce_{individual}-01"
    params:
        cutoff_umis = VALUES["cutoff_umis"],
        cutoff_doublets = VALUES["cutoff_doublets"]
    script:
        "scripts/01_remove_droplets.R"
        
# find outlier cells with low quality and remove them
rule remove_outliers:
    input: 
        sce_01 = rules.remove_droplets.output
    output:
        sce_02 = OUTPUT_DAT + "02_outl/{species}/sce_{individual}-02"
    params:
        cutoff_sum = VALUES["cutoff_sum"],
        cutoff_detected = VALUES["cutoff_detected"],
        cutoff_mitos = VALUES["cutoff_mitos"]
    script:
        "scripts/02_remove_outliers.R"
        
# normalize expression on sample level
rule normalize_expr:
    input: 
        sce_02 = rules.remove_outliers.output
    output:
        sce_03 = OUTPUT_DAT + "03_norm/{species}/sce_{individual}-03"
    script:
        "scripts/03_normalize_expr.R"
        
# extract hvgs and reduce dimensions for QC on sample level 
rule reduce_dims:
    input: 
        sce_03 = rules.normalize_expr.output
    output:
        sce_04 = OUTPUT_DAT + "04_dimr/{species}/sce_{individual}-04" 
    params:
        nr_hvgs = VALUES["nr_hvgs"]
    script:
        "scripts/04_reduce_dims.R"
        
#-------------------------------------------------------------------------------
# construct report files to monitor QC on sample level

rule make_sample_reports:
    input:
        sce_input = config["paths"]["cellranger_output"] + "/{species}/sce_{individual}-01",
        sce_01 = rules.remove_droplets.output,
        sce_02 = rules.remove_outliers.output,
        sce_03 = rules.normalize_expr.output,
        sce_04 = rules.reduce_dims.output
    output:
        OUTPUT_REP + "{species}/preprocessing_sample_report_{individual}.html"
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
        "preprocessing_sample_reports.Rmd" 
        
"""
Make one summary report on all files
Because this rule is dependend on previous outputs but this is not reflected
in input, it's generally turned off in config to avoid premature job submission
"""
print(config["run_preprocessing_summary"])
if config["run_preprocessing_summary"]:
  rule make_summary_report:
      input:
          sce_01_path = OUTPUT_DAT + "01_drop/"
      output:
          OUTPUT_REP + "preprocessing_summary.html"
      params:
          cutoff_sum = VALUES["cutoff_sum"],
          cutoff_detected = VALUES["cutoff_detected"],
          cutoff_mitos = VALUES["cutoff_mitos"],
          individuals = individuals,
          metadata = config["metadata"]["table"],
          color_tables = TABLES_PATH
      script:
          "preprocessing_summary.Rmd" 
