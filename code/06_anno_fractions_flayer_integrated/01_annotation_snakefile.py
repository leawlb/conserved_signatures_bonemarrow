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
species = get_list(metadata = METADATA, column = "Species_ID")
print(fractions)
clusters_hsc = list(range(1,14))
clusters_hsc = list(map(str,clusters_hsc))
clusters_hsc_useful = ["1", "2", "4", "5", "6", "7", "8", "9", "11", "12"]

clusters_str = list(range(1,18))
clusters_str = list(map(str,clusters_str))
clusters_str_useful = ["4", "6", "8", "9", "10", "13"]
# early mes, endos 1-3, fibros, CARs
# peris, smooth muscle and late CARs are too small clusters

# construct paths for all possible outputs/targets, required for rule all
targets = []

padj_cutoff = config["values"]["annotation"]["flayer_integrated_padj_cutoff"]
fc_cutoff = config["values"]["annotation"]["flayer_integrated_fc_cutoff"]

tf = "PC_" + str(padj_cutoff) + "_FC_" + str(fc_cutoff)
print(tf)

for c in clusters_hsc:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_hsc_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/hsc/markergenes/annotation_markergenes_hsc_cluster_" + c + ".html" ] 
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/hsc/quality_bulks/quality_bulks_hsc_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/hsc/quality_DGE/quality_dge_hsc_cluster_" + c + ".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/hsc/specific_DGE/" + tf + "/specific_hsc_cluster_" + c + ".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/hsc/shared_nDGE/" + tf + "/shared_hsc_cluster_" + c + ".html"]

for c in clusters_str:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_str_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/str/markergenes/annotation_markergenes_str_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/str/quality_bulks/quality_bulks_str_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/str/quality_DGE/quality_dge_str_cluster_" + c + ".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/str/specific_DGE/" + tf + "/specific_str_cluster_" + c + ".html"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/str/shared_nDGE/" + tf + "/shared_str_cluster_" + c + ".html"]

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
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/DGE/res_before_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/DGE/res_after_" + f + "-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/DGE/res_" + f + "_cluster_dfs-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/DGE/res_" + f + "_species_dfs-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/DGE/res_" + f + "_" + s + "_specific_df-11.csv" for s in species]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/nDGE/res_" + f + "_cluster_dfs-11"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/nDGE/res_" + f + "_species_dfs-11"]

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
        OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/{fraction}/markergenes/annotation_markergenes_{fraction}_cluster_{cluster}.html"
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
DESeq pre-processing and QC 
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
rule make_report_flayer_qc_bulks:
    input: 
        sce_11 = rules.assign_flayer_prelim.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        rld_11 = rules.qc_dge.output.rld_11,
        sva_11 = rules.qc_dge.output.sva_11,
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11
    output:
        OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/{fraction}/quality_bulks/quality_bulks_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_qc_dge.Rmd"

# get results from all pairwise comparisons
rule export_results_visualisation_before:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_before_11
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/DGE/res_before_{fraction}-11"
    script:
        "scripts/11_export_results_vis_flayer_integrated.R"    

# get results from pairwise comparisons for visualisation
rule export_results_visualisation_after:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_after_11
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/DGE/res_after_{fraction}-11"
    script:
        "scripts/11_export_results_vis_flayer_integrated.R"    

# progress report on DGE quality of DGE corrected for batch, age, SVs
rule make_report_flayer_pairwise_dge:
    input: 
        sce_11 = rules.assign_flayer_prelim.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11,
        res_before_11 =rules.export_results_visualisation_before.output,
        res_after_11 = rules.export_results_visualisation_after.output
    output:
        OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/{fraction}/quality_DGE/quality_dge_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_pairwise_dge.Rmd"
   
#-------------------------------------------------------------------------------
"""
Differential gene expression analysis 
"""  

# export shrunk FC results in convenient DF lists or csvs as required
rule export_results_dge:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_after_11 # after batch correction
    params:
        padj_cutoff = config["values"]["annotation"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["annotation"]["flayer_integrated_fc_cutoff"],
        species = species
    output:
        species_res_df = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/DGE/res_{fraction}_species_dfs-11",
        species_specific_csv = expand(OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/DGE/res_{{fraction}}_{species}_specific_df-11.csv", species = species),
        cluster_res_df = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/DGE/res_{fraction}_cluster_dfs-11"
    script:
        "scripts/11_export_results_dge_flayer_integrated.R"    

# report on species-specific genes
rule make_report_species_specific:
    input: 
        sce_11 = rules.assign_flayer_prelim.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq_11 = rules.preprocessing_dge.output.tdsq_after_11,
        cluster_res_df = rules.export_results_dge.output.cluster_res_df
    params:
        padj_cutoff = config["values"]["annotation"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["annotation"]["flayer_integrated_fc_cutoff"]
    output:
        OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/{fraction}/specific_DGE/" + tf + "/specific_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_specific.Rmd"

#-------------------------------------------------------------------------------
"""
Non-differential gene expression analysis (first layer)
To get shared genes between all species
"""  

# export results in convenient df lists or csvs
rule export_results_ndge:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_after_11
    params:
        padj_cutoff = config["values"]["annotation"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["annotation"]["flayer_integrated_fc_cutoff"],
        species = species
    output:
        species_res_df = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/nDGE/res_{fraction}_species_dfs-11",
        cluster_res_df = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/res_objects/" + tf + "/nDGE/res_{fraction}_cluster_dfs-11"
    script:
        "scripts/11_export_results_ndge_flayer_integrated.R"    
        
# report on shared genes for pval and cutoff selection
rule make_report_shared:
    input: 
        sce_11 = rules.assign_flayer_prelim.output,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq_11 = rules.preprocessing_dge.output.tdsq_after_11,
        cluster_res_df = rules.export_results_ndge.output.cluster_res_df
    params:
        padj_cutoff = config["values"]["annotation"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["annotation"]["flayer_integrated_fc_cutoff"]
    output:
        OUTPUT_BASE_PATH + "/reports/06_annotation/first_layer/{fraction}/shared_nDGE/" + tf + "/shared_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_shared.Rmd"

