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
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/02_drop/" + species + "/sce_" + ind + "-02"]
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/03_outl/" + species + "/sce_" + ind + "-03"]
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/04_norm/" + species + "/sce_" + ind + "-04"]
      targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/05_dimr/" + species + "/sce_" + ind + "-05"]
      targets = targets + [OUTPUT_BASE_PATH + "/reports/002_preprocessing/" + species + "/preprocessing_sample_report_" + ind + ".html"]
      
if config["run_preprocessing_summary"]:
  targets = targets + [OUTPUT_BASE_PATH + "/reports/002_preprocessing/preprocessing_summary.html"]

# local execution of non-demanding rules
localrules: all  

#-------------------------------------------------------------------------------

# define rules
rule all: # must contain all possible output paths from all rules below
    input:
        targets

# remove empty droplets and doublets to retain only single cells
rule remove_droplets:
    input: 
        sce_01 = OUTPUT_BASE_PATH + config["paths"]["cellranger_output"] + "/{species}/sce_{individual}-01"
    output:
        sce_02 = OUTPUT_BASE_PATH + "/sce_objects/02_drop/{species}/sce_{individual}-02"
    params:
        cutoff_umis = VALUES["cutoff_umis"],
        cutoff_doublets = VALUES["cutoff_doublets"]
    script:
        "scripts/02_remove_droplets.R"
        
# find outlier cells with low quality and remove them
rule remove_outliers:
    input: 
        sce_02 = rules.remove_droplets.output
    output:
        sce_03 = OUTPUT_BASE_PATH + "/sce_objects/03_outl/{species}/sce_{individual}-03"
    params:
        cutoff_sum = VALUES["cutoff_sum"],
        cutoff_detected = VALUES["cutoff_detected"],
        cutoff_mitos = VALUES["cutoff_mitos"]
    script:
        "scripts/03_remove_outliers.R"
        
# normalize expression on sample level
rule normalize_expr:
    input: 
        sce_03 = rules.remove_outliers.output
    output:
        sce_04 = OUTPUT_BASE_PATH + "/sce_objects/04_norm/{species}/sce_{individual}-04"
    script:
        "scripts/04_normalize_expr.R"
        
# extract hvgs and reduce dimensions for QC on sample level 
rule reduce_dims:
    input: 
        sce_04 = rules.normalize_expr.output
    output:
        sce_05 = OUTPUT_BASE_PATH + "/sce_objects/05_dimr/{species}/sce_{individual}-05" 
    params:
        nr_hvgs = VALUES["nr_hvgs"]
    script:
        "scripts/05_reduce_dims.R"
#-------------------------------------------------------------------------------
# construct report files to monitor QC on sample level

rule make_sample_reports:
    input:
        sce_01 = OUTPUT_BASE_PATH + config["paths"]["cellranger_output"] + "/{species}/sce_{individual}-01",
        sce_02 = rules.remove_droplets.output,
        sce_03 = rules.remove_outliers.output,
        sce_04 = rules.normalize_expr.output,
        sce_05 = rules.reduce_dims.output
    output:
        OUTPUT_BASE_PATH + "/reports/002_preprocessing/{species}/preprocessing_sample_report_{individual}.html"
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
        
# make one summary report on all files
"""
Because this rule is dependend on previous outputs but this is not reflected
in input:, it's generally turned off in config to avoid premature job submission
"""
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
