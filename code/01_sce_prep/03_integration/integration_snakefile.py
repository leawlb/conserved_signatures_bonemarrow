#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/01_sce_prep"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/01_sce_prep/03_integration"

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
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

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]

#-------------------------------------------------------------------------------

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/08_mrge/fractions/sce_" + f + "-08"]
  targets = targets + [OUTPUT_REP + "/fractions/integration_report_" + f + ".html"]
  for s in species:
    targets = targets + [OUTPUT_DAT + "/08_mrge/species/sce_" + s + "_" + f + "-08"]
    targets = targets + [OUTPUT_REP + "/species/integration_report_" + s + "_" + f + ".html"]

for s in species:
  for a in ages:
    targets = targets + [OUTPUT_DAT + "/08_mrge/ages/sce_" + s + "_" + a + "-08"]
    targets = targets + [OUTPUT_REP + "/ages/integration_report_" + s + "_" + a + ".html"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
# these scripts are almost identical with few exceptions.
# explanations in 08_merge_datasets_fraction.R

# Merge fractions (to be used for clustering, annotation)
rule merge_datasets_fractions:
    input: 
        sce_input_path = expand(OUTPUT_DAT + "/07_rfan/{s}/", s = species)
    output:
        sce_output_hsc = OUTPUT_DAT + "/08_mrge/fractions/sce_hsc-08",
        sce_output_str = OUTPUT_DAT + "/08_mrge/fractions/sce_str-08"
    params:
        individuals = individuals,
        species = species,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["nr_hvgs"],
        functions = "../../source/sce_functions.R"
    script:
        "scripts/08_merge_datasets_fractions.R" 
        
# Merge species 
rule merge_datasets_species:
    input: 
        sce_input_path = OUTPUT_DAT + "/07_rfan/{species}/"
    output:
        sce_output_hsc = OUTPUT_DAT + "/08_mrge/species/sce_{species}_hsc-08",
        sce_output_str = OUTPUT_DAT + "/08_mrge/species/sce_{species}_str-08"
    params:
        individuals = individuals,
        species = species,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["nr_hvgs"],
        functions = "../../source/sce_functions.R"
    script:
        "scripts/08_merge_datasets_species.R" 
        
# Merge ages (same as above but for each age separately) (for CCI)
rule merge_datasets_ages:
    input: 
        sce_input_path = OUTPUT_DAT + "/07_rfan/{species}/",
    output:
        sce_output_old = OUTPUT_DAT + "/08_mrge/ages/sce_{species}_old-08",
        sce_output_yng = OUTPUT_DAT + "/08_mrge/ages/sce_{species}_yng-08"
    params:
        individuals = individuals,
        species = species,
        ages = ages,
        samples_to_remove = config["samples_to_remove"],
        nr_hvgs = config["values"]["nr_hvgs"],
        functions = "../../source/sce_functions.R"
    script:
        "scripts/08_merge_datasets_ages.R" 
 
#-------------------------------------------------------------------------------
# Summary for fractions, species, and ages

rule make_reports_fractions:
    input: 
        sce_input = OUTPUT_DAT + "/08_mrge/fractions/sce_{fraction}-08"
    output:
        OUTPUT_REP + "/fractions/integration_report_{fraction}.html"
    params:
        colors_ref_path = COLORS_REF,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "integration_reports_fraction.Rmd" 
        
rule make_reports_species:
    input: 
        sce_input = OUTPUT_DAT + "/08_mrge/species/sce_{species}_{fraction}-08"
    output:
        OUTPUT_REP + "/species/integration_report_{species}_{fraction}.html"
    params:
        colors_ref_path = COLORS_REF,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "integration_reports_species.Rmd" 
        
rule make_reports_ages:
    input: 
        sce_input = OUTPUT_DAT + "/08_mrge/ages/sce_{species}_{age}-08"
    output:
        OUTPUT_REP + "/ages/integration_report_{species}_{age}.html"
    params:
        colors_ref_path = COLORS_REF,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "integration_reports_ages.Rmd" 
        
