#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_preparation"
OUTPUT_REP = OUTPUT_BASE + "/cci_objects/reports/01_cci_preparation"

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
for s in species:
  for a in age:
    targets = targets + [OUTPUT_DAT + "/10_ipis/ipi_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/11_idis/idi_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/11_idis/ilrs_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/11_idis/nrlrs_" + s + "_" + a]
    targets = targets + [OUTPUT_REP + "/cci_metrics_qc/cci_metrics_qc_report_" + s + "_" + a + ".html"]

targets = targets + [OUTPUT_REP +"/cci_metrics_ct_repertoire_summary.html"]
targets = targets + [OUTPUT_REP +"/cci_metrics_ct_interactome_summary.html"]
targets = targets + [OUTPUT_REP +"/cci_metrics_ctp_interactome_summary.html"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
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
        cci_input = OUTPUT_DAT + "/09_intl/interaction_list_{species}_{age}",
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05"
    output:
        ipi_list = OUTPUT_DAT + "/10_ipis/ipi_{species}_{age}"
    params:
        metrics_functions = "../../source/cci_functions_metrics.R"
    script:
        "scripts/10_get_ident_pair_info.R" 

# extract infos on identities
rule extract_idi:
    input:
        ipi_list = rules.extract_ipi.output
    output:
        idi_list = OUTPUT_DAT + "/11_idis/idi_{species}_{age}"
    params:
        metrics_functions = "../../source/cci_functions_metrics.R"
    script:
        "scripts/11_get_ident_info.R" 
        
# extract nr of ligands or receptors per identities
rule extract_nrlrs:
    input:
        cci_input = OUTPUT_DAT + "/09_intl/interaction_list_{species}_{age}",
        idi_list = rules.extract_idi.output
    output:
        ident_lrs_list = OUTPUT_DAT + "/11_idis/ilrs_{species}_{age}",
        ident_nrlrs_list = OUTPUT_DAT + "/11_idis/nrlrs_{species}_{age}"
    params:
        metrics_functions = "../../source/cci_functions_metrics.R"
    script:
        "scripts/11_get_ident_nrlrs_info.R" 
 
     
#-------------------------------------------------------------------------------

"""
Reports
"""

rule make_cci_qc_report_conditions:
    input:
        ident_pair_info = rules.extract_ipi.output.ipi_list,
        ident_info = rules.extract_idi.output.idi_list,
        ident_nrlrs_info = rules.extract_nrlrs.output.ident_nrlrs_list,
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05"
    output:
        OUTPUT_REP + "/cci_metrics_qc/cci_metrics_qc_report_{species}_{age}.html"
    params:
        colors_path = COLORS,
        sce_functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "cci_metrics_qc_report_conditions.Rmd" 

rule make_cci_ct_rep_sum:
    input:
        ident_nrlrs_info = expand(rules.extract_nrlrs.output.ident_nrlrs_list, species = species, age = age),
        sce_input = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05", species = species, age = age)
    output:
        OUTPUT_REP + "/cci_metrics_ct_repertoire_summary.html"
    params:
        colors_path = COLORS,
        sce_functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "cci_metrics_ct_repertoire_summary.Rmd" 
        
        
rule make_cci_ct_int_sum:
    input:
        ident_info = expand(rules.extract_idi.output.idi_list, species = species, age = age),
        sce_input = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05", species = species, age = age)
    output:
        OUTPUT_REP + "/cci_metrics_ct_interactome_summary.html"
    params:
        colors_path = COLORS,
        sce_functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "cci_metrics_ct_interactome_summary.Rmd" 
        
        
rule make_cci_ctp_int_sum:
    input:
        ident_pair_info = expand(rules.extract_ipi.output.ipi_list, species = species, age = age),
        sce_input = expand(OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05", species = species, age = age)
    output:
        OUTPUT_REP + "/cci_metrics_ctp_interactome_summary.html"
    params:
        colors_path = COLORS,
        sce_functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "cci_metrics_ctp_interactome_summary.Rmd" 
