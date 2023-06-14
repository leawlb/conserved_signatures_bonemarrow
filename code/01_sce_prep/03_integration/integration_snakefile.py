#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths and objects from config
OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])
TABLES_PATH = config["metadata"]["color_tables"]

# specific data and report output paths
OUTPUT_DAT = "/sce_objects/01_sce_prep/"
OUTPUT_REP = "/reports/01_sce_prep/03_integration/"

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

for f in fractions:
  targets = [OUTPUT_BASE + OUTPUT_DAT + "06_mrge/fractions/sce_" + f + "-06"]
  targets = targets + [OUTPUT_BASE + OUTPUT_DAT + "06_mrge/species/sce_" + s + "_" + f + "-06" for s in species]
  targets = targets + [OUTPUT_BASE + OUTPUT_REP + "species/integration_report_" + s + "_" + f + ".html" for s in species]

targets = targets + [OUTPUT_BASE + OUTPUT_REP + "fractions/integration_report_hsc.html"]
targets = targets + [OUTPUT_BASE + OUTPUT_REP + "fractions/integration_report_str.html"]
targets = targets + [OUTPUT_BASE + OUTPUT_REP + "species/integration_report_" + s + "_hsc.html" for s in species]
targets = targets + [OUTPUT_BASE + OUTPUT_REP + "species/integration_report_" + s + "_str.html" for s in species]

print(targets)
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
        
# Merge fractions 
rule merge_datasets_fractions:
    input: 
        sce_05_path = expand(OUTPUT_BASE + OUTPUT_DAT + "05_rfan/{s}/", s = species)
    output:
        sce_06_hsc = OUTPUT_BASE + OUTPUT_DAT + "06_mrge/fractions/sce_hsc-06",
        sce_06_str = OUTPUT_BASE + OUTPUT_DAT + "06_mrge/fractions/sce_str-06"
    params:
        individuals = individuals,
        species = species,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    script:
        "scripts/06_merge_datasets_fractions.R" 
        
# Merge species (same as above but for each species separately)
rule merge_datasets_species:
    input: 
        sce_05_path = OUTPUT_BASE + OUTPUT_DAT + "05_rfan/{species}/"
    output:
        sce_06_hsc = OUTPUT_BASE + OUTPUT_DAT + "06_mrge/species/sce_{species}_hsc-06",
        sce_06_str = OUTPUT_BASE + OUTPUT_DAT + "06_mrge/species/sce_{species}_str-06"
    params:
        individuals = individuals,
        species = species,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    script:
        "scripts/06_merge_datasets_species.R" 
 
# Summary for fractions
rule make_report_fractions:
    input: 
        sce_06 = OUTPUT_BASE + OUTPUT_DAT + "06_mrge/fractions/sce_{fraction}-06"
    output:
        OUTPUT_BASE + OUTPUT_REP + "fractions/integration_report_{fraction}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "integration_report_fraction.Rmd" 
        
# Summary for species
rule make_report_species:
    input: 
        sce_06 = OUTPUT_BASE + OUTPUT_DAT + "06_mrge/species/sce_{species}_{fraction}-06"
    output:
        OUTPUT_BASE + OUTPUT_REP + "species/integration_report_{species}_{fraction}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "integration_report_species.Rmd" 
        
