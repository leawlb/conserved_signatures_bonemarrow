#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT_01 = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own"
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/02_reclustering_other"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/02_reclustering_other"

RESOLUTION_OTHER = config["base"] + config["metadata_paths"]["resolution_other"]

RUN_PERM_GENESETS = config["run_permutation_genesets"]
RUN_PERM_BACKGROUND_SIGN = config["run_permutation_background_sign"]
RUN_PERM_BACKGROUND_MARK = config["run_permutation_background_mark"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
print(METADATA)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")

# define datasets to be tested
datasets_other_hsc = ["ts_hscs_progenitors", "ts_bone_marrow", "mus_weinreb_hspc", "mus_tm_bonemarrow", "nmr_sorted_hspc", "zeb_all_hspc"]
datasets_other_str = ["ts_all_stromal", "li_all_stromal", "mus_tik_stromal", "mus_bar_stromal"]
datasets_other = datasets_other_hsc + datasets_other_str
print(datasets_other)

#-------------------------------------------------------------------------------

targets = []
     
for d in datasets_other:
  targets = targets + [OUTPUT_DAT + "/03_recl/reclustered_" + d + "_list"]
  targets = targets + [OUTPUT_DAT + "/04_rcls/score_df_" + d + "_list"]
  targets = targets + [OUTPUT_REP + "/all/reclustering_other_report_" + d + ".html"]
  targets = targets + [OUTPUT_REP + "/final/reclustering_other_final_report_" + d + ".html"]
 
  # testing reclustering scores
  # targets = targets + [OUTPUT_REP + "/test_scores/test_reclustering_scores_" + d + ".html"]

  # permutation
  if RUN_PERM_GENESETS:
    targets = targets + [OUTPUT_DAT + "/05_perg/perm_score_df_mark_" + d]
    targets = targets + [OUTPUT_DAT + "/05_perg/perm_score_df_mmms_" + d]
    targets = targets + [OUTPUT_DAT + "/05_perg/perm_score_df_mmms_mark_" + d]
    targets = targets + [OUTPUT_REP + "/genesets/perm_genesets_" + d + ".html"]

  if RUN_PERM_BACKGROUND_SIGN:
    targets = targets + [OUTPUT_DAT + "/06_psig/perm_score_df_" + d]
    targets = targets + [OUTPUT_REP + "/conserved_signature/perm_conserved_signature_" + d + ".html"]

  if RUN_PERM_BACKGROUND_MARK:
    targets = targets + [OUTPUT_DAT + "/07_pmrk/perm_score_df_" + d]
    targets = targets + [OUTPUT_REP + "/conserved_markers/perm_conserved_markers_" + d + ".html"]


#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
      
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
       
"""
# RECLUSTERING
# 
# Reclustering datasets to test the ability of the conserved signatures
# to capture cell identity.
# Using our own datasets + from other species (human, zebrafish, naked mole rat)
# with:
#
# - conserved signatures
# - conserved marker genes
# - all BL6 marker genes 
# - random genes 
# - all conserved markers - random genes (other datasets only)
# - all BL6 markers - random genes (other datasets only)
#
# Additionally, permutation tests are performed.
"""

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING OTHER DATA
#
# Datasets from other species are prepared for analysis using 
# 00_prepare_datasets_snakefile.py
# All data generated from 00_prepare_datasets_snakefile.py is in metadata!
#
# Standard Seurat approach with standard options starting from raw counts
# but with aforementioned gene sets instead of HVGs.
"""
rule reclustering_other:
    input:
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/{dataset}",
        ensembl_sign = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_sign_{fraction}", fraction = fractions),
        ensembl_mark = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_mark_{fraction}", fraction = fractions),
        ensembl_mmms = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_mmms_{fraction}", fraction = fractions)
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"], 
        nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
        datasets_other_hsc = datasets_other_hsc,
        datasets_other_str = datasets_other_str
    output:
        seu_output = OUTPUT_DAT + "/03_recl/reclustered_{dataset}_list"
    script: 
        "02_scripts_other/03_reclustering_other.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING SCORES
# 
# In order to test how well a certain gene set was able to re-cluster a 
# dataset, we test different established (and new) metrics for
# comparing two clusterings, or for evaluating the quality of clustering.
# We use the re-clustered other datasets to test and evaluate each score.
# We choose to use following metrics for evaluating re-clustering:
#
# - adjusted rand index
# - variation of information
# - mean proportion of cells of a cell type per cluster 
#
# We use different scores because some of these scores increase or
# decrease with a higher clustering resolution/nr of clusters and we aim to 
# reduce this bias by using several different kinds of scores.
"""

rule test_reclustering_scores:
    input: 
        seu_list = rules.reclustering_other.output
    output:
        OUTPUT_REP + "/test_scores/test_reclustering_scores_{dataset}.html"
    conda:
        "../../envs/reclustering_scores.yml"
    script:
        "02_test_reclustering_scores.Rmd"


#-------------------------------------------------------------------------------  

"""
# REPORT OTHER DATASET SCORES
"""

# get the reclustering scores for re-clustered datasets from other species
rule reclustering_other_scores:
    input:
        seu_list = rules.reclustering_other.output,
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R"
    conda:
        "../../envs/reclustering_scores.yml"
    output:
        score_df_list = OUTPUT_DAT + "/04_rcls/score_df_{dataset}_list"
    script:
        "02_scripts_other/04_reclustering_other_scores.R"
        
# visualise re-clustering of other datasets, including with scores
rule reclustering_other_report:
    input:
        seu_list = rules.reclustering_other.output,
        score_df_list = rules.reclustering_other_scores.output
    output:
        OUTPUT_REP + "/all/reclustering_other_report_{dataset}.html"
    script: 
        "02_reclustering_other_report.Rmd"


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# PERMUTATION TESTS

# 1. Do other gene sets perform better than conserved signature genes?
#    Or is it just the larger number of genes?
# Permutation tests to compare:
# - conserved marker genes
# - all BL6 marker genes
# to gene sets comprised of conserved signature genes + x random genes 

# 2. Do conserved signature genes perform significantly better than random?
# Permutation tests to compare:
# - conserved signature genes
# to gene sets comprised of x random genes 

# Re-cluster seurat objects with random genes n = iteration times and compare
# the re-clustering scores of the original seurat object to the normal 
# distribution of the permuted scores.

# For permutation, the same resolution as the original gene set is chosen, 
# because each random set of genes might have a different optimal resolution. 

"""
#-------------------------------------------------------------------------------

if RUN_PERM_GENESETS:
  
  # requires conda channel genomedk
  rule permutation_genesets:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_sign = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_sign_{fraction}", fraction = fractions),
          ensembl_mark = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_mark_{fraction}", fraction = fractions),
          ensembl_mmms = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_mmms_{fraction}", fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          iterations = 20,
          nr_cores = 20,
          resolution_df = RESOLUTION_OTHER,
          #nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          #iterations = config["values"]["03_sce_analysis"]["iterations"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str     
      conda:
          "../../envs/reclustering_scores_permutation.yml"
      output:
          perm_score_df_mark = OUTPUT_DAT + "/05_perg/perm_score_df_mark_{dataset}",
          perm_score_df_mmms = OUTPUT_DAT + "/05_perg/perm_score_df_mmms_{dataset}",
          perm_score_df_mmms_mark = OUTPUT_DAT + "/05_perg/perm_score_df_mmms_mark_{dataset}"
      script: 
          "02_scripts_other/05_permutation_genesets.R"

  # visualise reclustering permutation
  rule permutation_genesets_report:
      input:
          score_df_list = rules.reclustering_other_scores.output,
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          perm_score_df_mark = rules.permutation_genesets.output.perm_score_df_mark,
          perm_score_df_mmms = rules.permutation_genesets.output.perm_score_df_mmms,
          perm_score_df_mmms_mark = rules.permutation_genesets.output.perm_score_df_mmms_mark,
          ensembl_sign = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_sign_{fraction}", fraction = fractions),
          ensembl_mark = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_mark_{fraction}", fraction = fractions),
          ensembl_mmms = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_mmms_{fraction}", fraction = fractions)
      params:
          resolution_df = RESOLUTION_OTHER,
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str
      output:
          OUTPUT_REP + "/genesets/perm_genesets_{dataset}.html"
      script: 
          "02_permutation_other_genesets_report.Rmd"
          
#-------------------------------------------------------------------------------

if RUN_PERM_BACKGROUND_SIGN:

  # requires conda channel genomedk
  rule permutation_background_sign:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_paths = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_sign_{fraction}", fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          iterations = 20,
          resolution_df = RESOLUTION_OTHER,
          #nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          #iterations = config["values"]["03_sce_analysis"]["iterations"],
          nr_cores = 20,
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          cons_level_use = "conserved_signature",
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str
      conda:
          "../../envs/reclustering_scores_permutation.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/06_psig/perm_score_df_{dataset}"
      script: 
          "02_scripts_other/06_permutation_background.R"


  # visualise reclustering permutation
  rule permutation_background_report_sign:
      input:
          seu_list = rules.reclustering_other.output,
          score_df_list = rules.reclustering_other_scores.output,
          perm_score_df = rules.permutation_background_sign.output
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "conserved_signature"
      output:
          OUTPUT_REP + "/conserved_signature/perm_conserved_signature_{dataset}.html"
      script: 
          "02_permutation_other_background_report.Rmd"

#-------------------------------------------------------------------------------

if RUN_PERM_BACKGROUND_MARK:

  # requires conda channel genomedk
  rule permutation_background_mark:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_paths = expand(OUTPUT_DAT_01 + "/02_endf/ensembl_mark_{fraction}", fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          iterations = 20,
          resolution_df = RESOLUTION_OTHER,
          #nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          #iterations = config["values"]["03_sce_analysis"]["nr_cores"],
          nr_cores = 20,
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          cons_level_use = "conserved_markers",
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str       
      conda:
          "../../envs/reclustering_scores_permutation.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/07_pmrk/perm_score_df_{dataset}"
      script: 
          "02_scripts_other/06_permutation_background.R"

  # visualise reclustering permutation
  rule permutation_background_report_mark:
      input:
          seu_list = rules.reclustering_other.output,
          score_df_list = rules.reclustering_other_scores.output,
          perm_score_df = rules.permutation_background_mark.output
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "conserved_markers"
      output:
          OUTPUT_REP + "/conserved_markers/perm_conserved_markers_{dataset}.html"
      script: 
          "02_permutation_other_background_report.Rmd"



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
FINAL REPORTS

# visualise one chosen re-clustering of other datasets with the best combination
# of scores
# for random features, always the same resolution as the conserved signature, 
# because each random set of genes (later for permutation) might have a dif-
# ferent optimal resolution. So choosing the same as the conserved signature
# is the easiest option

"""

rule reclustering_other_report_final:
    input:
        seu_list = rules.reclustering_other.output,
        score_df_list = rules.reclustering_other_scores.output,
        resolution_df = RESOLUTION_OTHER
    output:
        OUTPUT_REP + "/final/reclustering_other_final_report_{dataset}.html"
    script: 
      "02_reclustering_other_report_final.Rmd"
