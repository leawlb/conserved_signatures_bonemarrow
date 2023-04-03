#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
METADATA_PATH = config["metadata"]["table"]
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
IDENTIFIERS = config["metadata"]["identifiers"]
METADATA = pd.read_csv(config["metadata"]["table"])
VALUES =  config["values"]["preprocessing"]

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

Species_ID = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

# collect all paths for all possible outputs/targets, required for rule all
targets = []
for species in Species_ID:
  for ind in individuals:
    if species in ind:
      targets = targets + [OUTPUT_BASE_PATH + "/reports/003_check_mapping_general/" + species + "/mapping_sample_report_" + ind + ".html"]
      
#if config["run_preprocessing_summary"]:
#  targets = targets + [OUTPUT_BASE_PATH + "/reports/002_preprocessing/preprocessing_summary.html"]

# local execution of non-demanding rules
localrules: all  

#-------------------------------------------------------------------------------

# define rules
rule all: # must contain all possible output paths from all rules below
    input:
        targets

#-------------------------------------------------------------------------------
# construct report files to monitor QC on sample level

rule make_sample_reports:
    input:
        sce_03 = OUTPUT_BASE_PATH + "/sce_objects/03_outl/{species}/sce_{individual}-03",
        sce_fg = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/fourgenomes/sce_objects/01_cellranger_output/{species}/sce_{individual}-01"
    output:
        OUTPUT_BASE_PATH + "/reports/003_check_mapping_general/{species}/mapping_sample_report_{individual}.html"
    params:
        nr_hvgs = VALUES["nr_hvgs"],
        color_tables = TABLES_PATH
    script:
        # .Rmd files should be stored in snakefile working directory
        "mapping_sample_report.Rmd" 
        
# make one summary report on all files
"""
Because this rule is dependend on previous outputs but this is not reflected
in input:, it's generally turned off in config to avoid premature job submission

print(config["run_preprocessing_summary"])
if config["run_preprocessing_summary"]:
  rule make_summary_report:
      input:
          sce_02_path = OUTPUT_BASE_PATH + "/sce_objects/02_drop",
      output:
          OUTPUT_BASE_PATH + "/reports/002_preprocessing/preprocessing_summary.html"
      params:
          cutoff_sum = VALUES["cutoff_sum"],
          cutoff_detected = VALUES["cutoff_detected"],
          cutoff_mitos = VALUES["cutoff_mitos"],
          individuals = individuals,
          metadata = config["metadata"]["table"],
          color_tables = TABLES_PATH
      script:
          "preprocessing_summary.Rmd" 

"""
