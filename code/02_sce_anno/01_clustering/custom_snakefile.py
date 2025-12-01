#!/bin/python 

"""
Species-specific batch correction, clustering and custom annotation using 
louvain clustering. 
For comparison between species-specific and merged clustering 
to make sure that the merged clustering is not biased.
"""

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/02_sce_anno"

COLORS = config["base"] + config["metadata_paths"]["colors"]

VALUES =  config["values"]["02_sce_anno"]
BATCH_USE = VALUES["batch_use"] # which colData to use as batch

METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

#-------------------------------------------------------------------------------

targets = []

for f in fractions:
  targets = targets + [OUTPUT_REP + "/02_clustering/comparison_custom/custom_" + f + "_comparison_report.html"]
  targets = targets + [OUTPUT_REP + "/02_clustering/comparison_custom/custom_" + f + "_clustering_report.html"]
  
  for s in species:
    targets = targets + [OUTPUT_DAT + "/01_mnnc/comparison/sce_" + s + "_" + f + "_" + BATCH_USE + "-01"]
    targets = targets + [OUTPUT_REP + "/02_clustering/comparison_custom/report_" + s + "_" + f + ".html"]


#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# Batch correct and cluster species separately and compare to original
# clustering. Re-use the same scripts for a very basic, general 
# re-clustering using the same parameters.

# The batch correction has been performed in clustering_snakefile as a remainder
# from the basic species-specific clustering.
"""

#-------------------------------------------------------------------------------

"""
# Cluster species separately and compare to original clustering.
# Perform a custom, manual clustering to obtain a good
# cell type annotation for each species separately.
# Find the clustering parameters, then re-use the same scripts.

# Do this in an .Rmd file so the same time I check all the plots
# and add the annotation I can also save the object immediately.
# This is not good practice but it's the most convenient option

# One custom .Rdm per condition (species + fraction)
"""

# MMUS

rule louvain_clustering_mmus_hsc:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mmus_hsc_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mmus_hsc.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mmus_hsc.Rmd"

rule louvain_clustering_mmus_str:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mmus_str_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mmus_str.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mmus_str.Rmd"

# MCAS

rule louvain_clustering_mcas_hsc:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mcas_hsc_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mcas_hsc.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mcas_hsc.Rmd"

rule louvain_clustering_mcas_str:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mcas_str_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mcas_str.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mcas_str.Rmd"

# MSPR

rule louvain_clustering_mspr_hsc:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mspr_hsc_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mspr_hsc.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mspr_hsc.Rmd"

rule louvain_clustering_mspr_str:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mspr_str_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mspr_str.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mspr_str.Rmd"

# MCAR

rule louvain_clustering_mcar_hsc:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mcar_hsc_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mcar_hsc.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mcar_hsc.Rmd"

rule louvain_clustering_mcar_str:
    input: 
        sce_species = OUTPUT_DAT + "/01_mnnc/comparison/sce_mcar_str_" + BATCH_USE + "-01"
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/report_mcar_str.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "custom_clustering/custom_mcar_str.Rmd"


"""
# REPORTS comparing between species, and between the original merged 
# and the custom species-specific clustering.
"""

# HSC reports

# TODO: add input later
rule custom_hsc_comparison_report:
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/custom_hsc_comparison_report.html"
    resources:
        mem_mb = 70000,
        queue = "medium-debian"
    threads: 10
    script:
        "custom_clustering/custom_hsc_comparison_report.Rmd"

# TODO: add input later
rule custom_hsc_clustering_report:
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/custom_hsc_clustering_report.html"
    resources:
        mem_mb = 70000,
        queue = "medium-debian"
    conda:
        "../../envs/ggalluvial.yml"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    threads: 10
    script:
        "custom_clustering/custom_hsc_clustering_report.Rmd"

# STR reports

# TODO: add input later
rule custom_str_comparison_report:
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/custom_str_comparison_report.html"
    resources:
        mem_mb = 70000,
        queue = "medium-debian"
    threads: 10
    script:
        "custom_clustering/custom_str_comparison_report.Rmd"

# TODO: add input later
rule custom_str_clustering_report:
    output:
        OUTPUT_REP + "/02_clustering/comparison_custom/custom_str_clustering_report.html"
    resources:
        mem_mb = 70000,
        queue = "medium-debian"
    conda:
        "../../envs/ggalluvial.yml"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    threads: 10
    script:
        "custom_clustering/custom_str_clustering_report.Rmd"

