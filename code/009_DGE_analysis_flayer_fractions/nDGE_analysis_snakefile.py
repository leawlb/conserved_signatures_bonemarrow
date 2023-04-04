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

clusters_hsc = ["1", "2", "4", "5", "6", "7", "8", "9", "11", "12"]
clusters_str = ["4", "6", "8", "9", "10", "13"]

# construct paths for all possible outputs/targets, required for rule all
targets = []

padj_cutoff = config["values"]["DGE"]["flayer_integrated_padj_cutoff"]
fc_cutoff = config["values"]["DGE"]["flayer_integrated_fc_cutoff"]
tf = "PC_" + str(padj_cutoff) + "_FC_" + str(fc_cutoff)
print(tf)

for c in clusters_hsc:
  targets = targets + [OUTPUT_BASE_PATH + "/reports/009_DGE_analysis_flayer/hsc/shared_nDGE/" + tf + "/shared_hsc_cluster_" + c + ".html"]

for c in clusters_str:
  targets = targets + [OUTPUT_BASE_PATH + "/reports/009_DGE_analysis_flayer/str/shared_nDGE/" + tf + "/shared_str_cluster_" + c + ".html"]

for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/nDGE/res_" + f + "_cluster_dfs-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/nDGE/res_" + f + "_species_dfs-11"]

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
Non-differential gene expression analysis (first layer)
To get shared genes between all species
"""  

# export results in convenient df lists or csvs
rule export_results_ndge:
    input:
        tdsq = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_after_{fraction}-11" # after batch correction
    params:
        padj_cutoff = config["values"]["DGE"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["DGE"]["flayer_integrated_fc_cutoff"],
        species = species
    output:
        species_res_df = OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/nDGE/res_{fraction}_species_dfs-11",
        cluster_res_df = OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/" + tf + "/nDGE/res_{fraction}_cluster_dfs-11"
    script:
        "scripts/13_export_results_ndge_flayer_integrated.R"    
        
# report on shared genes for pval and cutoff selection
rule make_report_shared:
    input: 
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_after_{fraction}-11", # after batch correction
        cluster_res_df = rules.export_results_ndge.output.cluster_res_df,
        gene_list = config["metadata"]["gene_list_hsc"],
    params:
        padj_cutoff = config["values"]["DGE"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["DGE"]["flayer_integrated_fc_cutoff"]
    output:
        OUTPUT_BASE_PATH + "/reports/009_DGE_analysis_flayer/{fraction}/shared_nDGE/" + tf + "/shared_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_shared.Rmd"


