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
  targets = targets + [OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/hsc/quality_bulks/quality_bulks_hsc_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/hsc/quality_DGE/quality_dge_hsc_cluster_" + c + ".html"]

for c in clusters_str:
  targets = targets + [OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/str/quality_bulks/quality_bulks_str_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/str/quality_DGE/quality_dge_str_cluster_" + c + ".html"]

for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/dsq_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/rld_objects/results_rld/rld_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/sva_objects/results_sva/sva_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_before_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_after_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/res_objects/DGE/res_before_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/res_objects/DGE/res_after_" + f + "-11"]

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
DESeq pre-processing and QC 
"""

# form pseudobulks/aggregate per cluster per sample
rule aggregate_convert:
    input:
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
    output:
        dsq_11 = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/dsq_{fraction}-11"
    script:
        "scripts/12_aggregate_flayer_integrated.R"
    
# produce rlog-transformed values, identify hidden sources of variations (SV) for visualisation
rule qc_dge:
    input:
        dsq_11 = rules.aggregate_convert.output
    output:
        rld_11 = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/rld_objects/results_rld/rld_{fraction}-11",
        sva_11 = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/sva_objects/results_sva/sva_{fraction}-11"
    script:
        "scripts/12_qc_dge_flayer_integrated.R"  
    
# preprocessing and transformation for DGE
rule preprocessing_dge:
    input:
        dsq_11 = rules.aggregate_convert.output,
        sva_11 = rules.qc_dge.output.sva_11
    output:
        tdsq_before_11 = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_before_{fraction}-11",
        tdsq_after_11 = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_after_{fraction}-11"
    script:
        "scripts/12_prepro_dge_flayer_integrated.R"    

# progress report on pseudo-bulk quality, SVs and effect of BC with SVs
rule make_report_flayer_qc_bulks:
    input: 
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        rld_11 = rules.qc_dge.output.rld_11,
        sva_11 = rules.qc_dge.output.sva_11,
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11
    output:
        OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/{fraction}/quality_bulks/quality_bulks_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_qc_dge.Rmd"

# get results from all pairwise comparisons
rule export_results_visualisation_before:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_before_11
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/res_objects/DGE/res_before_{fraction}-11"
    script:
        "scripts/12_export_results_vis_flayer_integrated.R"    

# get results from pairwise comparisons for visualisation
rule export_results_visualisation_after:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_after_11
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/res_objects/DGE/res_after_{fraction}-11"
    script:
        "scripts/12_export_results_vis_flayer_integrated.R"    

# progress report on DGE quality of DGE corrected for batch, age, SVs
rule make_report_flayer_pairwise_dge:
    input: 
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11,
        res_before_11 =rules.export_results_visualisation_before.output,
        res_after_11 = rules.export_results_visualisation_after.output
    output:
        OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/{fraction}/quality_DGE/quality_dge_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_pairwise_dge.Rmd"
   
