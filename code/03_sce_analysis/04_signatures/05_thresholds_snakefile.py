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
targets = targets + [OUTPUT_DAT + "/ts_bone_marrow_recl_by_p.rds"]
targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_recl_by_t.rds"]
targets = targets + [OUTPUT_DAT + "/ts_bone_marrow_recl_by_t.rds"]

targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_scores_by_p.rds"]
targets = targets + [OUTPUT_DAT + "/ts_bone_marrow_scores_by_p.rds"]
targets = targets + [OUTPUT_DAT + "/ts_hscs_progenitors_scores_by_t.rds"]
targets = targets + [OUTPUT_DAT + "/ts_bone_marrow_scores_by_t.rds"]

targets = targets + [OUTPUT_REP + "/05_threshold_reclustering_report_ts_hscs_progenitors_p.html"]
targets = targets + [OUTPUT_REP + "/05_threshold_reclustering_report_ts_bone_marrow_p.html"]
targets = targets + [OUTPUT_REP + "/05_threshold_reclustering_report_ts_hscs_progenitors_t.html"]
targets = targets + [OUTPUT_REP + "/05_threshold_reclustering_report_ts_bone_marrow_t.html"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

# recluster ts_hscs_progenitors using different sets of signatureso obtained
# after sliding the threshold on the nDGE analysis p value (by_p)

# I select ts_hspcs_progenitors first, because it looks relatively stable in
# figure S4 even with different numbers of genes (no big jumps or bumps
# of scores with changes in resolution, because it is one of the first 
# human datasets and because the annotation has a better resolution than 
# the ts adult whole bone marrow

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# reclustering after sliding p value threshold

# ts_hscs_progenitors
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

# ts_bone_marrow
rule recl_ts_bone_marrow_by_p:
    input: 
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_bone_marrow",
        ensembl_hum = ENSEMBL_HUM,
        ensembl_mus = ENSEMBL_MUS,
        metadata_path = METADATA_PATH
    output:
        seu_recl_output = OUTPUT_DAT + "/ts_bone_marrow_recl_by_p.rds"
    params:
        threshold_path = OUTPUT_DAT,
        dataset = "ts_bone_marrow",
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        nr_cores = 10
    resources:
        mem_mb = 80000,
        queue = "medium-debian"
    threads: 10
    script:
        "05_thresholds/0X_recluster_by_p.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# reclustering after sliding findmarkers log2FC

# ts_hscs_progenitors
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

# ts_bone_marrow
rule recl_ts_bone_marrow_by_t:
    input: 
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_bone_marrow",
        ensembl_hum = ENSEMBL_HUM,
        ensembl_mus = ENSEMBL_MUS,
        metadata_path = METADATA_PATH
    output:
        seu_recl_output = OUTPUT_DAT + "/ts_bone_marrow_recl_by_t.rds"
    params:
        threshold_path = OUTPUT_DAT,
        dataset = "ts_bone_marrow",
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        nr_cores = 10
    resources:
        mem_mb = 80000,
        queue = "medium-debian"
    threads: 10
    script:
        "05_thresholds/0X_recluster_by_t.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# get the reclustering scores 
# based on 02_scripts_other/04_reclustering_other_scores.R
# after sliding nDGE p values

# ts_hscs_progenitors
rule reclustering_other_scores_ts_hscs_progenitors_p:
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

# ts_bone_marrow
rule reclustering_other_scores_ts_bone_marrow_p:
    input:
        seu_list = rules.recl_ts_bone_marrow_by_p.output,
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        threshold_used = "conserved_signature_treshold_testing_by_p"
    conda:
        "../../envs/reclust_scores_perm_others.yml" 
    output:
        score_df_list = OUTPUT_DAT + "/ts_bone_marrow_scores_by_p.rds"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_thresholds/0X_recl_scores.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# get the reclustering scores 
# based on 02_scripts_other/04_reclustering_other_scores.R
# after sliding log2FC FindMarkers
# using the same script

# ts_hscs_progenitors
rule reclustering_other_scores_ts_hscs_progenitors_t:
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

# ts_bone_marrow
rule reclustering_other_scores_ts_bone_marrow_t:
    input:
        seu_list = rules.recl_ts_bone_marrow_by_t.output,
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        threshold_used = "conserved_signature_treshold_testing_by_t"
    conda:
        "../../envs/reclust_scores_perm_others.yml" 
    output:
        score_df_list = OUTPUT_DAT + "/ts_bone_marrow_scores_by_t.rds"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_thresholds/0X_recl_scores.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# report after sliding ndge p values

# ts_hscs_progenitors
rule threshold_reclustering_report_ts_hscs_progenitors_p:
    input:
        score_df_list = rules.reclustering_other_scores_ts_hscs_progenitors_p.output,
    params:
        threshold_used = "conserved_signature_treshold_testing_by_p",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        dataset = "ts_hscs_progenitors"
    output:
        OUTPUT_REP + "/05_threshold_reclustering_report_ts_hscs_progenitors_p.html"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_threshold_reclustering_report.Rmd"

# ts_bone_marrow
rule threshold_reclustering_report_ts_bone_marrow_p:
    input:
        score_df_list = rules.reclustering_other_scores_ts_bone_marrow_p.output,
    params:
        threshold_used = "conserved_signature_treshold_testing_by_p",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        dataset = "ts_bone_marrow"
    output:
        OUTPUT_REP + "/05_threshold_reclustering_report_ts_bone_marrow_p.html"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_threshold_reclustering_report.Rmd"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# report after sliding findmarkers log2fc

# ts_hscs_progenitors
rule threshold_reclustering_report_ts_hscs_progenitors_t:
    input:
        score_df_list = rules.reclustering_other_scores_ts_hscs_progenitors_t.output,
    params:
        threshold_used = "conserved_signature_treshold_testing_by_t",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        dataset = "ts_hscs_progenitors"
    output:
        OUTPUT_REP + "/05_threshold_reclustering_report_ts_hscs_progenitors_t.html"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_threshold_reclustering_report.Rmd"

# ts_bone_marrow
rule threshold_reclustering_report_ts_bone_marrow_t:
    input:
        score_df_list = rules.reclustering_other_scores_ts_bone_marrow_t.output,
    params:
        threshold_used = "conserved_signature_treshold_testing_by_t",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        dataset = "ts_bone_marrow"
    output:
        OUTPUT_REP + "/05_threshold_reclustering_report_ts_bone_marrow_t.html"
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 4
    script:
        "05_threshold_reclustering_report.Rmd"