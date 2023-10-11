#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/03_cci_analysis"
OUTPUT_REP = OUTPUT_BASE + "/reports/06_cci_analysis/02_cci_ed"

COLORS = config["base"] + config["metadata_paths"]["colors"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")

#VALUES = config["values"]["05_cci_analysis"]

#-------------------------------------------------------------------------------

targets = []
for a in age:
  targets = targets + [OUTPUT_DAT + "/02_ed/ages/cci_dist_" + a]
#  targets = targets + [OUTPUT_REP + "/ages/pca_report_" + a + ".html"]
#  targets = targets + [OUTPUT_REP + "/ages/pca_explore_" + a + ".html"]

for s in species:
  targets = targets + [OUTPUT_DAT + "/02_ed/species/cci_dist_" + s]
  targets = targets + [OUTPUT_REP + "/species/report_ed_" + s + ".html"]
#  targets = targets + [OUTPUT_REP + "/species/pca_explore_" + s + ".html"]
  
targets = targets + [OUTPUT_DAT + "/02_ed/all/cci_dist_all"]
targets = targets + [OUTPUT_REP + "/all/report_ed_all.html"]
#targets = targets + [OUTPUT_REP + "/all/pca_explore_all.html"]

#-------------------------------------------------------------------------------


localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

rule calculate_ed_ages:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/ages/interaction_list_{age}"
    output:
        cci_dist_output = OUTPUT_DAT + "/02_ed/ages/cci_dist_{age}"
    script:
        "scripts/01_ED.R"  
 
rule calculate_ed_species:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/species/interaction_list_{species}"
    output:
        cci_dist_output = OUTPUT_DAT + "/02_ed/species/cci_dist_{species}"
    script:
        "scripts/01_ED.R"         

rule calculate_ed_all:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/all/interaction_list_all"
    output:
        cci_dist_output = OUTPUT_DAT + "/02_ed/all/cci_dist_all"
    script:
        "scripts/01_ED.R"    

#-------------------------------------------------------------------------------       
 
rule report_ed_species:
    input: 
        cci_dist_input = rules.calculate_ed_species.output,
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/species/interaction_list_{species}"
    output:
        OUTPUT_REP + "/species/report_ed_{species}.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
    script:
        "cci_ed_report.Rmd"         

rule report_ed_all:
    input: 
        cci_dist_input = rules.calculate_ed_all.output,
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/all/interaction_list_all"
    output:
        OUTPUT_REP + "/all/report_ed_all.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
    script:
        "cci_ed_report.Rmd"         
