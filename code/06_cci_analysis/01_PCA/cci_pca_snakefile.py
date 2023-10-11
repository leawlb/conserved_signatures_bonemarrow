#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/03_cci_analysis"
OUTPUT_REP = OUTPUT_BASE + "/reports/06_cci_analysis/01_cci_pca"

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
  targets = targets + [OUTPUT_DAT + "/01_pca/ages/pca_list_" + a]
  targets = targets + [OUTPUT_REP + "/ages/pca_report_" + a + ".html"]
  targets = targets + [OUTPUT_REP + "/ages/pca_explore_" + a + ".html"]

for s in species:
  targets = targets + [OUTPUT_DAT + "/01_pca/species/pca_list_" + s]
  targets = targets + [OUTPUT_REP + "/species/pca_report_" + s + ".html"]
  targets = targets + [OUTPUT_REP + "/species/pca_explore_" + s + ".html"]
  
targets = targets + [OUTPUT_DAT + "/01_pca/all/pca_list_all"]
targets = targets + [OUTPUT_REP + "/all/pca_report_all.html"]
targets = targets + [OUTPUT_REP + "/all/pca_explore_all.html"]

#-------------------------------------------------------------------------------


localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

rule calculate_pca_ages:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/ages/interaction_list_{age}"
    output:
        pca_output = OUTPUT_DAT + "/01_pca/ages/pca_list_{age}"
    params:
        cci_functions_analysis = "../../source/cci_functions_analysis.R"
    script:
        "scripts/01_calculate_pca.R"  
 
rule calculate_pca_species:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/species/interaction_list_{species}"
    output:
        pca_output = OUTPUT_DAT + "/01_pca/species/pca_list_{species}"
    params:
        cci_functions_analysis = "../../source/cci_functions_analysis.R"
    script:
        "scripts/01_calculate_pca.R"         

rule calculate_pca_all:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/all/interaction_list_all"
    output:
        pca_output = OUTPUT_DAT + "/01_pca/all/pca_list_all"
    params:
        cci_functions_analysis = "../../source/cci_functions_analysis.R"
    script:
        "scripts/01_calculate_pca.R"    
        
#-------------------------------------------------------------------------------

rule report_pca_ages:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/ages/interaction_list_{age}",
        pca_input = rules.calculate_pca_ages.output,
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/ages/pca_report_{age}.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    script:
        "cci_pca_report.Rmd"  

rule report_pca_species:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/species/interaction_list_{species}",
        pca_input = rules.calculate_pca_species.output,
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/species/pca_report_{species}.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    script:
        "cci_pca_report.Rmd"  
        
rule report_pca_all:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/all/interaction_list_all",
        pca_input = rules.calculate_pca_all.output,
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/all/pca_report_all.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    script:
        "cci_pca_report.Rmd"  

#-------------------------------------------------------------------------------
     

rule explore_pca_ages:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/ages/interaction_list_{age}",
        pca_input = rules.calculate_pca_ages.output,
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/ages/pca_explore_{age}.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    script:
        "cci_pca_explore_report.Rmd"  

rule explore_pca_species:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/species/interaction_list_{species}",
        pca_input = rules.calculate_pca_species.output,
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/species/pca_explore_{species}.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    script:
        "cci_pca_explore_report.Rmd"  
        
rule explore_pca_all:
    input: 
        cci_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/12_mrge/all/interaction_list_all",
        pca_input = rules.calculate_pca_all.output,
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/all/pca_explore_all.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    script:
        "cci_pca_explore_report.Rmd"  

