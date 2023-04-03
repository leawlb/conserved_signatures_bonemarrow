#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import sys 

# paths from config
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["table"])
print(OUTPUT_BASE_PATH)

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")
print(fractions)
clusters_hsc = list(range(1,14))
clusters_hsc = list(map(str,clusters_hsc))

clusters_str = list(range(1,18))
clusters_str = list(map(str,clusters_str))

padj_cutoff = config["values"]["DGE"]["flayer_integrated_padj_cutoff"]
fc_cutoff = config["values"]["DGE"]["flayer_integrated_fc_cutoff"]
tf = "PC_" + str(padj_cutoff) + "_FC_" + str(fc_cutoff)
print(tf)

# construct paths for all possible outputs/targets, required for rule all
targets = []

for c in clusters_hsc:
  targets = targets + [OUTPUT_BASE_PATH + "/reports/009_DGE_analysis_flayer/hsc/specific_DGE/" + tf + "/specific_hsc_cluster_" + c + ".html"]

for c in clusters_str:
  targets = targets + [OUTPUT_BASE_PATH + "/reports/009_DGE_analysis_flayer/str/specific_DGE/" + tf + "/specific_str_cluster_" + c + ".html"]

for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/DGE/res_" + f + "_cluster_dfs-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/DGE/res_" + f + "_species_dfs-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/DGE/res_" + f + "_" + s + "_specific_df-11.csv" for s in species]

#print(targets)

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
    cluster="[0-9]+",

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------
"""
Differential gene expression analysis 
"""  

# export shrunk FC results in convenient DF lists or csvs as required
rule export_results_dge:
    input:
        tdsq = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_after_{fraction}-11" # after batch correction
    params:
        padj_cutoff = config["values"]["DGE"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["DGE"]["flayer_integrated_fc_cutoff"],
        species = species
    output:
        species_res_df = OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/DGE/res_{fraction}_species_dfs-11",
        species_specific_csv = expand(OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/DGE/res_{{fraction}}_{species}_specific_df-11.csv", species = species),
        cluster_res_df = OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/DGE/res_{fraction}_cluster_dfs-11"
    script:
        "scripts/13_export_results_dge_flayer_integrated.R"    

# report on species-specific genes
rule make_report_species_specific:
    input: 
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq_11 = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_after_{fraction}-11",
        cluster_res_df = rules.export_results_dge.output.cluster_res_df
    params:
        padj_cutoff = config["values"]["DGE"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["DGE"]["flayer_integrated_fc_cutoff"]
    output:
        OUTPUT_BASE_PATH + "/reports/009_DGE_analysis_flayer/{fraction}/specific_DGE/" + tf + "/specific_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_specific.Rmd"
