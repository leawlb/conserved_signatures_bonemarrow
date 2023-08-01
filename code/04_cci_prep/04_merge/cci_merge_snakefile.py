#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

# paths and objects from config
OUTPUT_BASE = config["base"] + config["data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_preparation"

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
ages = get_list(metadata = METADATA, column = "Age_ID")

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]

#-------------------------------------------------------------------------------

targets = []

for a in ages:
  targets = targets + [OUTPUT_DAT + "/12_mrge/ages/interaction_list_" + a]
  
for s in species:
  targets = targets + [OUTPUT_DAT + "/12_mrge/species/interaction_list_" + s]

targets = targets + [OUTPUT_DAT + "/12_mrge/all/interaction_list_all"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------

# Merge to ages
rule merge_ccio_ages:
    input: 
        cci_input_path = expand(OUTPUT_DAT + "/09_intl/interaction_list_{s}_{{age}}/", s = species)
    output:
        cci_output = OUTPUT_DAT + "/12_mrge/ages/interaction_list_{age}"
    script:
        "scripts/12_merge_ccio.R" 
   
# merge to species     
rule merge_ccio_species:
    input: 
        cci_input_path = expand(OUTPUT_DAT + "/09_intl/interaction_list_{{species}}_{a}/", a = ages)
    output:
        cci_output = OUTPUT_DAT + "/12_mrge/species/interaction_list_{species}"
    script:
        "scripts/12_merge_ccio.R" 
 
# merge all cci objects for total comparison
rule merge_ccio_all:
    input: 
        cci_input_path = expand(OUTPUT_DAT + "/09_intl/interaction_list_{s}_{a}/", a = ages, s = species)
    output:
        cci_output = OUTPUT_DAT + "/12_mrge/all/interaction_list_all"
    script:
        "scripts/12_merge_ccio.R"        
