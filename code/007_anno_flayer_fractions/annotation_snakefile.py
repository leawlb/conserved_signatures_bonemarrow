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

# construct paths for all possible outputs/targets, required for rule all
targets = []

for c in clusters_hsc:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_hsc_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/007_annotation/hsc/markergenes/annotation_markergenes_hsc_cluster_" + c + ".html" ] 

for c in clusters_str:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_str_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/007_annotation/str/markergenes/annotation_markergenes_str_cluster_" + c + ".html"] 

for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_markers/markers_" + f]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_go/go_" + f]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_" + f + "-11"]

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
    cluster="[0-9]+",

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------
# dummy rule to allow using clusters as wildcard
output = expand(OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_hsc_cluster_{cluster}-sep", cluster = clusters_hsc)
output = output + expand(OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_str_cluster_{cluster}-sep", cluster = clusters_str)
rule separate_sce:
    input: 
        sce_10 = expand(OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10", fraction = fractions)
    output:
        sce_11_sep = output
    script:
        "scripts/00_separate_dummy.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
annotation for first layer of broad clusters (flayer)
"""

# get marker genes for each cluster and do go analysis on the marker genes.
rule find_markers:
    input: 
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:
        results_markers = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_markers/markers_{fraction}",
    params:
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
    script:
        "scripts/11_markers_flayer.R"

# perform preliminary GO analysis for overview (Very simple)
rule go_analysis:
    input: 
        results_markers = rules.find_markers.output
    output:
        results_go = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/results_go/go_{fraction}"
    script:
        "scripts/11_go_flayer.R"

# progress report on marker gene expression and GO 
rule make_report_flayer_markergenes:
    input: 
        results_markers = rules.find_markers.output,
        results_go = rules.go_analysis.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:
        OUTPUT_BASE_PATH + "/reports/007_annotation/{fraction}/markergenes/annotation_markergenes_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_markergenes.Rmd"
        
# assign labels 
rule assign_flayer:
    input:
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:      
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11"
    script:
        "scripts/11_anno_flayer.R"
      
