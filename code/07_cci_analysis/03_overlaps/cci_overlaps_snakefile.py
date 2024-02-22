#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/03_cci_analysis"
OUTPUT_REP = OUTPUT_BASE + "/cci_objects/reports/03_cci_analysis/03_overlaps"

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

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_REP + "/ages/overlap_ages_ct_interactome.html"]
targets = targets + [OUTPUT_REP + "/ages/overlap_ages_ctp_interactome.html"]
targets = targets + [OUTPUT_REP + "/ages/overlap_ages_ct_repertoire.html"]
#targets = targets + [OUTPUT_REP + "/ages/overlap_species_nrlrs_report.html"]

print(targets)

#-------------------------------------------------------------------------------


localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------


rule overlap_ages_ct_interactome:
    input: 
        cci_met_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/11_idis/idi_{species}_{age}", species = species, age = age)
    output:
        OUTPUT_REP + "/ages/overlap_ages_ct_interactome.html"
    params: 
        species = species,
        age = age,
        colors_path = COLORS,
        colors = "../../source/colors.R",
    script:
        "overlap_ages_ct_interactome.Rmd"  

rule overlap_ages_ctp_interactome:
    input: 
        cci_met_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/10_ipis/ipi_{species}_{age}", species = species, age = age)
    output:
        OUTPUT_REP + "/ages/overlap_ages_ctp_interactome.html"
    params: 
        species = species,
        age = age,
        colors_path = COLORS,
        colors = "../../source/colors.R",
    script:
        "overlap_ages_ctp_interactome.Rmd"  


rule overlap_ages_ct_repertoire:
    input: 
        cci_met_ilrs_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/11_idis/ilrs_{species}_{age}", species = species, age = age),
        cci_met_idi_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/11_idis/idi_{species}_{age}", species = species, age = age),
        sce_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05", species = species, age = age),
        perc_df_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/06_perc/perc_df_{species}_{age}", species = species, age = age),
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/ages/overlap_ages_ct_repertoire.html"
    params: 
        species = species,
        age = age,
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    script:
        "overlap_ages_ct_repertoire.Rmd"  
"""

rule overlap_species_nrlrs_report:
    input: 
        cci_met_ilrs_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/11_idis/ilrs_{species}_yng", species = species),
        sce_paths = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_yng-05", species = species),
        sce_input_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_input_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10",
    output:
        OUTPUT_REP + "/ages/overlap_species_nrlrs_report.html"
    params: 
        species = species,
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R",
    script:
        "overlap_species_nrlrs.Rmd" 

"""  
