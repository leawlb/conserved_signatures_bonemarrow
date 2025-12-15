#!/bin/python 

"""
TODO LATER
"""

#-------------------------------------------------------------------------------

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
METADATA_PATH = config["base"] + config["metadata_paths"]["metadata"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/05_thresholds"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/05_thresholds"

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"]
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"]

COLORS = config["base"] + config["metadata_paths"]["colors"]

#-------------------------------------------------------------------------------

targets = []
targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_p.rds"]
targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_t.rds"]

targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_scores_by_p.rds"]
targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_scores_by_t.rds"]

targets = targets + [OUTPUT_REP + "/05_threshold_reclustering_report_p.html"]
targets = targets + [OUTPUT_REP + "/05_threshold_reclustering_report_t.html"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

# recluster ts_hscs_progenitors using different sets of signatureso obtained
# after sliding the threshold on the nDGE analysis p value (by_p)

# I select ts_hspcs_progenitors, because it looks relatively stable in
# figure S4 even with different numbers of genes (no big jumps or bumps
# of scores with changes in resolution, because it is one of the first 
# human datasets and because the annotation has a better resolution than 
# the ts adult whole bone marrow


rule recl_ts_hscs_progenitors_by_p:
    input: 
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_hscs_progenitors",
        ensembl_hum = ENSEMBL_HUM,
        ensembl_mus = ENSEMBL_MUS,
        metadata_path = METADATA_PATH
    output:
        seu_recl_output = OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_p.rds"
    params:
        threshold_path = OUTPUT_DAT,
        dataset = "ts_hscs_progenitors",
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        nr_cores = 10
    resources:
        mem_mb = 80000,
        queue = "medium-debian"
    threads: 10
    script:
        "05_thresholds/0X_recluster_by_p.R"

rule recl_ts_hscs_progenitors_by_t:
    input: 
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_hscs_progenitors",
        ensembl_hum = ENSEMBL_HUM,
        ensembl_mus = ENSEMBL_MUS,
        metadata_path = METADATA_PATH
    output:
        seu_recl_output = OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_t.rds"
    params:
        threshold_path = OUTPUT_DAT,
        dataset = "ts_hscs_progenitors",
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        nr_cores = 10
    resources:
        mem_mb = 80000,
        queue = "medium-debian"
    threads: 10
    script:
        "05_thresholds/0X_recluster_by_t.R"

# get the reclustering scores 
# based on 02_scripts_other/04_reclustering_other_scores.R
rule reclustering_other_scores_p:
    input:
        seu_list = rules.recl_ts_hscs_progenitors_by_p.output,
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        threshold_used = "conserved_signature_treshold_testing_by_p"
    conda:
        "../../envs/reclust_scores_perm_others.yml" 
    output:
        score_df_list = OUTPUT_DAT + "/ts_hscs_progenitors_scores_by_p.rds"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_thresholds/0X_recl_scores.R"

# using the same script
rule reclustering_other_scores_t:
    input:
        seu_list = rules.recl_ts_hscs_progenitors_by_t.output,
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        threshold_used = "conserved_signature_treshold_testing_by_t"
    conda:
        "../../envs/reclust_scores_perm_others.yml" 
    output:
        score_df_list = OUTPUT_DAT + "/ts_hscs_progenitors_scores_by_t.rds"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_thresholds/0X_recl_scores.R"

rule threshold_reclustering_report_p:
    input:
        score_df_list = rules.reclustering_other_scores_p.output,
    params:
        threshold_used = "conserved_signature_treshold_testing_by_p",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        dataset = "ts_hscs_progenitors"
    output:
        OUTPUT_REP + "/05_threshold_reclustering_report_p.html"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_threshold_reclustering_report.Rmd"

rule threshold_reclustering_report_t:
    input:
        score_df_list = rules.reclustering_other_scores_t.output,
    params:
        threshold_used = "conserved_signature_treshold_testing_by_t",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        dataset = "ts_hscs_progenitors"
    output:
        OUTPUT_REP + "/05_threshold_reclustering_report_t.html"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_threshold_reclustering_report.Rmd"