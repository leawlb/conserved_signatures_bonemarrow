#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import sys 

# paths from config
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["raw"])
print(OUTPUT_BASE_PATH)

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")
print(fractions)
clusters_hsc = list(range(1,14))
clusters_hsc = list(map(str,clusters_hsc))

clusters_str = list(range(1,18))
clusters_str = list(map(str,clusters_str))

# construct paths for all possible outputs/targets, required for rule all
targets = []

targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_hsc_cluster_" + c + "-sep" for c in clusters_hsc]
targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_str_cluster_" + c + "-sep" for c in clusters_str]
targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/hsc/annotation_flayer_markergenes_hsc_cluster_" + c + ".html" for c in clusters_hsc] 
targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/str/annotation_flayer_markergenes_str_cluster_" + c + ".html" for c in clusters_str] 
targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/hsc/annotation_flayer_quality_dge_hsc_cluster_" + c + ".html" for c in clusters_hsc] 
targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/str/annotation_flayer_quality_dge_str_cluster_" + c + ".html" for c in clusters_str] 
targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/hsc/annotation_flayer_pairwise_dge_hsc_cluster_" + c + ".html" for c in clusters_hsc]
targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/str/annotation_flayer_pairwise_dge_str_cluster_" + c + ".html" for c in clusters_str]

for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_markers/markers_" + f]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_go/go_" + f]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_go/go_" + f]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_prelimanno_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/dsq_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_rld/rld_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_sva/sva_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/tdsq_before_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/tdsq_after_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/res_before_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/res_after_" + f + "-11"]

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
FIRST LAYER, INTEGRATED DATASET
"""

#-------------------------------------------------------------------------------
"""
preliminary annotation for first layer of broad clusters (flayer)
"""

# get marker genes for each cluster and do go analysis on the marker genes.
rule find_markers:
    input: 
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:
        results_markers = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_markers/markers_{fraction}",
    params:
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
    script:
        "scripts/11_markers_flayer_integrated.R"

# perform preliminary GO analysis for overview
rule go_analysis:
    input: 
        results_markers = rules.find_markers.output
    output:
        results_go = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_go/go_{fraction}"
    script:
        "scripts/11_go_flayer_integrated.R"

# progress report on marker gene expression and GO for first layer, integrated
rule make_report_flayer_markergenes:
    input: 
        results_markers = rules.find_markers.output,
        results_go = rules.go_analysis.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:
        OUTPUT_BASE_PATH + "/reports/06_annotation/{fraction}/annotation_flayer_markergenes_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_markergenes.Rmd"
        
# assign preliminary labels for overview 
rule assign_flayer_prelim:
    input:
        sce_10 = OUTPUT_BASE_PATH + "/sce_objects/10_lcls_fractions/sce_{fraction}-10"
    output:      
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_prelimanno_{fraction}-11"
    script:
        "scripts/11_prelimanno_flayer_integrated.R"
      
#-------------------------------------------------------------------------------
"""
Differential gene expression (DGE) analysis of first layer
"""

# form pseudobulks/aggregate per cluster per sample
rule aggregate_convert:
    input:
        sce_11 = rules.assign_flayer_prelim.output
    output:
        dsq_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/dsq_{fraction}-11"
    script:
        "scripts/11_aggregate_flayer_integrated.R"
    
# produce rlog-transformed values, identify hidden sources of variations (SV) for visualisation
rule qc_dge:
    input:
        dsq_11 = rules.aggregate_convert.output
    output:
        rld_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_rld/rld_{fraction}-11",
        sva_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/other/results_sva/sva_{fraction}-11"
    script:
        "scripts/11_qc_dge_flayer_integrated.R"  
    
# preprocessing and transformation for DGE
rule preprocessing_dge:
    input:
        dsq_11 = rules.aggregate_convert.output,
        sva_11 = rules.qc_dge.output.sva_11
    output:
        tdsq_before_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/tdsq_before_{fraction}-11",
        tdsq_after_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/tdsq_after_{fraction}-11"
    script:
        "scripts/11_prepro_dge_flayer_integrated.R"    

# progress report on pseudo-bulk quality, SVs and effect of BC with SVs
rule make_report_flayer_qc_dge:
    input: 
        sce_11 = rules.assign_flayer_prelim.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        rld_11 = rules.qc_dge.output.rld_11,
        sva_11 = rules.qc_dge.output.sva_11,
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11
    output:
        OUTPUT_BASE_PATH + "/reports/06_annotation/{fraction}/annotation_flayer_quality_dge_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_qc_dge.Rmd"

# get results from all pairwise comparisons
rule pairwise_comparisons_dge_before:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_before_11
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/res_before_{fraction}-11"
    script:
        "scripts/11_pairwisecomp_dge_flayer_integrated.R"    

# get results from all pairwise comparisons
rule pairwise_comparisons_dge_after:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_after_11
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/dsq_objects/res_after_{fraction}-11"
    script:
        "scripts/11_pairwisecomp_dge_flayer_integrated.R"    

# progress report on pairwise comparisons between species, corrected for batch, age, SVs
rule make_report_flayer_pairwise_dge:
    input: 
        sce_11 = rules.assign_flayer_prelim.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11,
        res_before_11 =rules.pairwise_comparisons_dge_before.output,
        res_after_11 = rules.pairwise_comparisons_dge_after.output
    output:
        OUTPUT_BASE_PATH + "/reports/06_annotation/{fraction}/annotation_flayer_pairwise_dge_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_pairwise_dge.Rmd"
   
    
    
    
    
