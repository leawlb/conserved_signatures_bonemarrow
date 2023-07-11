#!/bin/python 

# Prepping SCEs for CCI calculation

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_preparation/"
OUTPUT_REP = OUTPUT_BASE + "/reports/04_cci_calc/"

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

#targets = targets + [OUTPUT_REP + "sce_downsampling_report.html"]
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
        interaction_list = OUTPUT_DAT + "08_intl/interaction_list_{species}_{age}"
    output:
        ident_pair_info = OUTPUT_DAT + "09_ipis/ipi_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/10_get_ident_pair_info.R" 

# extract infos on identities
rule extract_idi:
    input:
        ident_pair_info = rules.extract_ipi.output
    output:
        ident_info = OUTPUT_DAT + "10_idis/idi_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/11_get_ident_info.R" 
        
# extract nr of ligands or receptors per identities
rule extract_nrlrs:
    input:
        interaction_list = OUTPUT_BASE + "/cci_objects/02_cci_calc/04_intl/interaction_list_{species}_{age}",
        ident_pair_info = rules.extract_idi.output
    output:
        ident_lrs_info = OUTPUT_DAT + "10_idis/ilrs_{species}_{age}",
        ident_nrlrs_info = OUTPUT_DAT + "10_idis/nrlrs_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/11_get_ident_nrlrs_info.R" 
 
     
#-------------------------------------------------------------------------------

"""
Reports
"""


