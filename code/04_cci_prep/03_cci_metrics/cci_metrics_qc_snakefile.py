#!/bin/python 

# Prepping SCEs for CCI calculation

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])
TABLES_PATH = config["metadata"]["color_tables"]

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_preparation/"
OUTPUT_REP = OUTPUT_BASE + "/reports/04_cci_prep/01_cci_preparation/"

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")

targets = []

for s in species:
  for a in age:
    targets = targets + [OUTPUT_DAT + "10_ipis/ipi_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "11_idis/idi_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "11_idis/ilrs_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "11_idis/nrlrs_" + s + "_" + a]
    targets = targets + [OUTPUT_REP + "cci_metrics_qc/cci_metrics_qc_report_" + s + "_" + a + ".html"]

targets = targets + [OUTPUT_REP + "cci_metrics_qc/cci_metrics_qc_summary.html"]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
 
#-------------------------------------------------------------------------------

"""
Metrics

These functions are in principle similar to PerCellQC metrics.
For each identity pair and each identity, important metrics such as 
number of interaction and number of ligands or receptors are extracted and
stored separately

"""

# extract infos on identity pairs
rule extract_ipi:
    input:
        interaction_list = OUTPUT_DAT + "09_intl/interaction_list_{species}_{age}"
    output:
        ident_pair_info = OUTPUT_DAT + "10_ipis/ipi_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/10_get_ident_pair_info.R" 

# extract infos on identities
rule extract_idi:
    input:
        ident_pair_info = rules.extract_ipi.output
    output:
        ident_info = OUTPUT_DAT + "11_idis/idi_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/11_get_ident_info.R" 
        
# extract nr of ligands or receptors per identities
rule extract_nrlrs:
    input:
        interaction_list = OUTPUT_DAT + "09_intl/interaction_list_{species}_{age}",
        ident_pair_info = rules.extract_idi.output
    output:
        ident_lrs_info = OUTPUT_DAT + "11_idis/ilrs_{species}_{age}",
        ident_nrlrs_info = OUTPUT_DAT + "11_idis/nrlrs_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/11_get_ident_nrlrs_info.R" 
 
     
#-------------------------------------------------------------------------------

"""
Reports
"""

rule make_cci_qc_report_conditions:
    input:
        ident_pair_info = rules.extract_ipi.output.ident_pair_info,
        ident_info = rules.extract_idi.output.ident_info,
        ident_nrlrs_info = rules.extract_nrlrs.output.ident_nrlrs_info,
        sce_input = OUTPUT_DAT + "05_down/sce_{species}_{age}-05"
    output:
        OUTPUT_REP + "cci_metrics_qc/cci_metrics_qc_report_{species}_{age}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "cci_metrics_qc_report_conditions.Rmd" 


rule make_cci_qc_summary:
    input:
        ident_pair_info = expand(rules.extract_ipi.output.ident_pair_info, species = species, age = age),
        ident_info = expand(rules.extract_idi.output.ident_info, species = species, age = age),
        ident_nrlrs_info = expand(rules.extract_nrlrs.output.ident_nrlrs_info, species = species, age = age),
        sce_input = expand(OUTPUT_DAT + "05_down/sce_{species}_{age}-05", species = species, age = age)
    output:
        OUTPUT_REP + "cci_metrics_qc/cci_metrics_qc_summary.html"
    params:
        color_tables = TABLES_PATH
    script:
        "cci_metrics_qc_summary.Rmd" 
        
