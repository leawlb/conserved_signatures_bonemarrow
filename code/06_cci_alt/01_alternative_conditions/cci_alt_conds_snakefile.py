#!/bin/python 

# Monster Snakefile of entire CCI process to compare conditions.
# For more details, check out 04_cci_prep and 06_cci_analysis.

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/02_cci_alt/01_conditions"
OUTPUT_REP = OUTPUT_BASE + "/cci_objects/reports/02_cci_alt/01_conditions"

COLORS = config["base"] + config["metadata_paths"]["colors"]

LRDB_OUT =  config["base"] + config["metadata_paths"]["lrdb_out"]
ASSIGNMENT = config["base"] + config["metadata_paths"]["assignment"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")
conds = ["main", "nctf", "nlvl", "ndwn", "nprp"]

VALUES = config["values"]["05_cci_prep"]


#-------------------------------------------------------------------------------

targets = []
for s in species:
  for a in age:
    for c in conds:
      targets = targets + [OUTPUT_DAT + "/01_calc_cond/cci_" + s + "_" + a + "_" + c]
      targets = targets + [OUTPUT_DAT + "/02_mtrc_cond/cci_metrics_" + s + "_" + a + "_" + c]
      targets = targets + [OUTPUT_DAT + "/03_cpca_cond/cci_pca_" + s + "_" + a + "_" + c]

targets = targets + [OUTPUT_REP + "/conditions_summary.html"]
    
#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------


"""

Calculate CCIs under different conditions:
  
 - all pre-processing steps (main)
 - no cutoff
 - no rank levelling
 - no downsampling
 - no pre-processing steps
 
The following rules summarize the entire process of CCI calculation.

"""

rule calculate_main:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05",
        lrdb = LRDB_OUT
    output:
        cci_output = OUTPUT_DAT + "/01_calc_cond/cci_{species}_{age}_main"
    params:
        main_functions = "../../source/cci_functions_calculation.R",
        alternative_functions = "../../source/cci_functions_alt_conds.R",
        min_perc = VALUES["min_perc"],
        top_level = VALUES["top_level"],
        assay_use = "downsampled",
        condition = "all preprocessing steps"
    script:
        "scripts/01_calculate_conds.R" 

rule calculate_nocutoff:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05",
        lrdb = LRDB_OUT
    output:
        cci_output = OUTPUT_DAT + "/01_calc_cond/cci_{species}_{age}_nctf"
    params:
        main_functions = "../../source/cci_functions_calculation.R",
        alternative_functions = "../../source/cci_functions_alt_conds.R",
        min_perc = 0,
        top_level = VALUES["top_level"],
        assay_use = "downsampled",
        condition = "no cutoff"
    script:
        "scripts/01_calculate_conds.R" 
        
rule calculate_nolevel:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05",
        lrdb = LRDB_OUT
    output:
        cci_output = OUTPUT_DAT + "/01_calc_cond/cci_{species}_{age}_nlvl"
    params:
        main_functions = "../../source/cci_functions_calculation.R",
        alternative_functions = "../../source/cci_functions_alt_conds.R",
        min_perc = VALUES["min_perc"],
        top_level = "0",
        assay_use = "downsampled",
        condition = "no rank level"
    script:
        "scripts/01_calculate_conds.R" 

rule calculate_nodownsampling:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05",
        lrdb = LRDB_OUT
    output:
        cci_output = OUTPUT_DAT + "/01_calc_cond/cci_{species}_{age}_ndwn"
    params:
        main_functions = "../../source/cci_functions_calculation.R",
        alternative_functions = "../../source/cci_functions_alt_conds.R",
        min_perc = VALUES["min_perc"],
        top_level = VALUES["top_level"],
        assay_use = "counts",
        condition = "no downsampling"
    script:
        "scripts/01_calculate_conds.R" 
        
rule calculate_nopreprocessing:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05",
        lrdb = LRDB_OUT
    output:
        cci_output = OUTPUT_DAT + "/01_calc_cond/cci_{species}_{age}_nprp"
    params:
        main_functions = "../../source/cci_functions_calculation.R",
        alternative_functions = "../../source/cci_functions_alt_conds.R",
        min_perc = 0,
        top_level = 0,
        assay_use = "counts",
        condition = "no pre-processing steps"
    script:
        "scripts/01_calculate_conds.R" 


"""

Calculate CCI metrics under different conditions:
  
The following rules summarizes the entire process of CCI metrics calculation.

"""

rule metrics:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05",
        cci_input =  OUTPUT_DAT + "/01_calc_cond/cci_{species}_{age}_{condition}"
    output:
        metrics_lists = OUTPUT_DAT + "/02_mtrc_cond/cci_metrics_{species}_{age}_{condition}"
    params:
        metrics_functions = "../../source/cci_functions_metrics.R"
    script:
        "scripts/02_metrics_conds.R" 

"""

Do PCA under different conditions:
  
"""

rule pca:
    input:
        cci_input =  OUTPUT_DAT + "/01_calc_cond/cci_{species}_{age}_{condition}"
    output:
        pca_output = OUTPUT_DAT + "/03_cpca_cond/cci_pca_{species}_{age}_{condition}"
    params:
        metrics_functions = "../../source/cci_functions_metrics.R"
    script:
        "scripts/03_pca_conds.R"
        
        
"""

Reports

"""

rule cond_report:
    input:
        metrics_path = expand(rules.metrics.output.metrics_lists, age = age, species = species, condition = conds),
        pca_path = expand(rules.pca.output.pca_output, age = age, species = species, condition = conds)
    output:
        OUTPUT_REP + "/conditions_summary.html"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"    
    script:
        "cci_alt_conds_report.Rmd" 
