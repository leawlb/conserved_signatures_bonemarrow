#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths and objects from config
OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])
TABLES_PATH = config["metadata"]["color_tables"]

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/01_sce_prep/"
OUTPUT_REP = OUTPUT_BASE + "/reports/01_sce_prep/04_integration/"

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
ages = get_list(metadata = METADATA, column = "Age_ID")

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "08_mrge/fractions/sce_" + f + "-08"]
  targets = targets + [OUTPUT_REP + "fractions/integration_report_" + f + ".html"]
  
for f in fractions:
  for s in species:
    targets = targets + [OUTPUT_DAT + "08_mrge/species/sce_" + s + "_" + f + "-08"]
    targets = targets + [OUTPUT_REP + "species/integration_report_" + s + "_" + f + ".html"]

for s in species:
  for a in ages:
    targets = targets + [OUTPUT_DAT + "08_mrge/ages/sce_" + s + "_" + a + "-08"]
    targets = targets + [OUTPUT_REP + "ages/integration_report_" + s + "_" + a + ".html"]

#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
        
# Merge fractions (for annotation)
rule merge_datasets_fractions:
    input: 
        sce_05_path = expand(OUTPUT_DAT + "05_rfan/{s}/", s = species),
        dmg_list = OUTPUT_DAT + "07_dmgs/dmgs_list_all"
    output:
        sce_08_hsc = OUTPUT_DAT + "08_mrge/fractions/sce_hsc-08",
        sce_08_str = OUTPUT_DAT + "08_mrge/fractions/sce_str-08"
    params:
        individuals = individuals,
        species = species,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    script:
        "scripts/08_merge_datasets_fractions.R" 
        
# Merge species (same as above but for each species separately) (for NMF)
rule merge_datasets_species:
    input: 
        sce_05_path = OUTPUT_DAT + "05_rfan/{species}/",
        dmg_list = OUTPUT_DAT + "07_dmgs/dmgs_list_all"
    output:
        sce_08_hsc = OUTPUT_DAT + "08_mrge/species/sce_{species}_hsc-08",
        sce_08_str = OUTPUT_DAT + "08_mrge/species/sce_{species}_str-08"
    params:
        individuals = individuals,
        species = species,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    script:
        "scripts/08_merge_datasets_species.R" 
        
# Merge ages (same as above but for each age separately) (for CCI)
rule merge_datasets_ages:
    input: 
        sce_05_path = OUTPUT_DAT + "05_rfan/{species}/",
        dmg_list = OUTPUT_DAT + "07_dmgs/dmgs_list_all"
    output:
        sce_08_old = OUTPUT_DAT + "08_mrge/ages/sce_{species}_old-08",
        sce_08_yng = OUTPUT_DAT + "08_mrge/ages/sce_{species}_yng-08"
    params:
        individuals = individuals,
        species = species,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        ages = ages
    script:
        "scripts/08_merge_datasets_ages.R" 
 
#-------------------------------------------------------------------------------

# Summary for fractions
rule make_report_fractions:
    input: 
        sce_08 = OUTPUT_DAT + "08_mrge/fractions/sce_{fraction}-08"
    output:
        OUTPUT_REP + "fractions/integration_report_{fraction}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "integration_report_fraction.Rmd" 
        
# Summary for species
rule make_report_species:
    input: 
        sce_08 = OUTPUT_DAT + "08_mrge/species/sce_{species}_{fraction}-08"
    output:
        OUTPUT_REP + "species/integration_report_{species}_{fraction}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "integration_report_species.Rmd" 
        
# Summary for ages
rule make_report_ages:
    input: 
        sce_08 = OUTPUT_DAT + "08_mrge/ages/sce_{species}_{age}-08"
    output:
        OUTPUT_REP + "ages/integration_report_{species}_{age}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "integration_report_ages.Rmd" 
        
