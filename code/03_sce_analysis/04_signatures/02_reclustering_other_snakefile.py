#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT_OWN = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own"
OUTPUT_DAT_OTHER = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/02_reclustering_other"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/02_reclustering_other"

RESOLUTION_OTHER = config["base_input"] + config["metadata_paths"]["resolution_other"]
FRAME_OTHER = config["base_input"] + config["metadata_paths"]["frame_other"]

# permutation: each gene set individually vs random 
RUN_SIGN_RAND_OTHER_PERM = config["run_sign_rand_other_permutation"]
RUN_MARK_RAND_OTHER_PERM = config["run_mark_rand_other_permutation"] # set to FALSE
RUN_MMMS_RAND_OTHER_PERM = config["run_mmms_rand_other_permutation"]

# permutation: marker sets vs signature + random
RUN_GNST_SIGN_RAND_OTHER_PERM = config["run_genesets_sign_rand_other_permutation"]
print(RUN_GNST_SIGN_RAND_OTHER_PERM)

# run pval correction AT THE END
RUN_PVAL_CORRECTION = config["run_pval_correction"]

METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")

# define datasets to be tested
# nmr_sorted_hspc, zeb_all_hspc
datasets_other_hsc = ["ts_hscs_progenitors", "ts_bone_marrow", "mus_weinreb_hspc", "mus_tm_bonemarrow"]
datasets_other_str = ["ts_all_stromal", "li_all_stromal", "mus_tik_stromal", "mus_bar_stromal"]
datasets_other = datasets_other_hsc + datasets_other_str
print(datasets_other)

#-------------------------------------------------------------------------------

targets = []
     
for d in datasets_other:
  targets = targets + [OUTPUT_DAT_OTHER + "/03_recl/reclustered_" + d + "_list"]
  targets = targets + [OUTPUT_DAT_OTHER + "/04_rcls/score_df_" + d + "_list"]
  targets = targets + [OUTPUT_REP + "/all/reclustering_other_report_" + d + ".html"]
  targets = targets + [OUTPUT_REP + "/final/reclustering_other_final_report_" + d + ".html"]
 
  #testing reclustering scores
  #targets = targets + [OUTPUT_REP + "/test_scores/test_reclustering_scores_" + d + ".html"]

  if RUN_SIGN_RAND_OTHER_PERM:
    targets = targets + [OUTPUT_DAT_OTHER + "/05_psig/perm_score_df_" + d]
    targets = targets + [OUTPUT_REP + "/signature/perm_sign_rand_" + d + ".html"]
    targets = targets + [OUTPUT_DAT_OTHER + "/08_expp/sign-vs-rand_" + d]

  if RUN_MARK_RAND_OTHER_PERM:
    targets = targets + [OUTPUT_DAT_OTHER + "/06_pmrk/perm_score_df_" + d]
    targets = targets + [OUTPUT_REP + "/conserved_markers/perm_cons_rand_" + d + ".html"]
    targets = targets + [OUTPUT_DAT_OTHER + "/08_expp/mark-vs-rand_" + d]

  if RUN_MMMS_RAND_OTHER_PERM:
    targets = targets + [OUTPUT_DAT_OTHER + "/07_pmms/perm_score_df_" + d]
    targets = targets + [OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-rand_" + d]


  # permutation
  if RUN_GNST_SIGN_RAND_OTHER_PERM:
    targets = targets + [OUTPUT_DAT_OTHER + "/07_perg/perm_score_df_mark_" + d]
    targets = targets + [OUTPUT_DAT_OTHER + "/07_perg/perm_score_df_mmms_" + d]
    targets = targets + [OUTPUT_DAT_OTHER + "/07_perg/perm_score_df_mmms_mark_" + d]
    targets = targets + [OUTPUT_REP + "/genesets/perm_genesets_" + d + ".html"]
    targets = targets + [OUTPUT_DAT_OTHER + "/08_expp/mark-vs-signrand_" + d]
    targets = targets + [OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-signrand_" + d]
    #targets = targets + [OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-markrand_" + d]

if RUN_PVAL_CORRECTION:
  targets = targets + [OUTPUT_DAT_OTHER + "/09_crpv/all_corrected_pval"]

#-------------------------------------------------------------------------------

localrules: all, pval_correction, sign_vs_rand, mark_vs_signrand, mmms_vs_signrand, mmms_vs_rand, mark_vs_rand

rule all: 
    input:
        targets
      
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
       
"""
# RECLUSTERING
# 
# Reclustering datasets to test the ability of the signature
# to capture cell identity.
# Using our own datasets + from other species (human, zebrafish, naked mole rat)
# with:
#
# - signatures
# - conserved marker genes
# - all BL6 marker genes 
# - random genes 
# - all conserved markers - random genes (other datasets only) # still true?
# - all BL6 markers - random genes (other datasets only) # still true?
#
# Additionally, permutation tests are performed.
"""

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

print(config["base"] + config["metadata_paths"]["datasets_other_path"])

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
        ensembl_sign = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_sign_{fraction}", fraction = fractions),
        ensembl_mark = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mark_{fraction}", fraction = fractions),
        ensembl_mmms = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mmms_{fraction}", fraction = fractions)
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"], 
        nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
        datasets_other_hsc = datasets_other_hsc,
        datasets_other_str = datasets_other_str
    output:
        seu_output = OUTPUT_DAT_OTHER + "/03_recl/reclustered_{dataset}_list"
    resources:
          mem_mb=50000,
          queue="long"
    threads: 4
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
# - mean proportion of cells of a cell type per cluster (top 10% ranked)
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
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
    conda:
        "../../envs/reclust_scores_perm_others.yml"
    resources:
          mem_mb=50000,
          queue="long"
    threads: 4
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
        "../../envs/reclust_scores_perm_others.yml" 
    output:
        score_df_list = OUTPUT_DAT_OTHER + "/04_rcls/score_df_{dataset}_list"
    resources:
          mem_mb=50000,
          queue="medium"
    threads: 4
    script:
        "02_scripts_other/04_reclustering_other_scores.R"
        
# visualise re-clustering of other datasets, including with scores
rule reclustering_other_report:
    input:
        seu_list = rules.reclustering_other.output,
        score_df_list = rules.reclustering_other_scores.output
    output:
        OUTPUT_REP + "/all/reclustering_other_report_{dataset}.html"
    resources:
          mem_mb=50000,
          queue="medium"
    threads: 4
    script: 
        "02_reclustering_other_report.Rmd"


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# PERMUTATION TESTS


# 1. Do signature genes perform significantly better than random?
# Permutation tests to compare:
# - signature genes
# - conserved markers
# - BL6 markers
# to gene sets comprised of x random genes 

# 2. Do other gene sets perform better than signature genes?
#    Or is it just the larger number of genes?
# Permutation tests to compare:
# - conserved marker genes
# - BL6 marker genes
# to gene sets comprised of signature genes + x random genes 


# Re-cluster seurat objects with random genes n = iteration times and compare
# the re-clustering scores of the original seurat object to the normal 
# distribution of the permuted scores.

# For permutation, the same resolution as the original gene set is chosen, 
# because each random set of genes might have a different optimal resolution. 

"""
#-------------------------------------------------------------------------------

if RUN_SIGN_RAND_OTHER_PERM:

  rule perm_sign_rand_other:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_paths = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_sign_{fraction}", fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          resolution_df = RESOLUTION_OTHER,
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          iterations = config["values"]["03_sce_analysis"]["iterations"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          cons_level_use = "conserved_signature", # termed only "signature" now
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str
      conda:
          "../../envs/reclust_scores_perm_others.yml"
      output:
          perm_score_df = OUTPUT_DAT_OTHER + "/05_psig/perm_score_df_{dataset}"
      resources:
          mem_mb=200000,
          queue="long"
      threads: config["values"]["03_sce_analysis"]["nr_cores"]
      script: 
          "02_scripts_other/05_permutation_other_rand.R"

  # export pvals: compare signature to random background
  rule sign_vs_rand:
      input:
          orig_score_df_list_input = rules.reclustering_other_scores.output.score_df_list,
          perm_score_df_input = rules.perm_sign_rand_other.output.perm_score_df
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "conserved_signature",
          cons_level_use_alt = "seu_sign", 
          comparison = "sign-vs-rand"
      output:
          pval_score_df_output = OUTPUT_DAT_OTHER + "/08_expp/sign-vs-rand_{dataset}"
      resources:
          mem_mb=20000,
          queue="medium"
      threads: 4
      script:
          "02_scripts_other/08_export_pval.R"
          
  # visualise reclustering permutation
  rule permutation_report_sign_rand:
      input:
          seu_list = rules.reclustering_other.output,
          score_df_list = rules.reclustering_other_scores.output,
          perm_score_df = rules.perm_sign_rand_other.output
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "conserved_signature"
      output:
          OUTPUT_REP + "/signature/perm_sign_rand_{dataset}.html"
      resources:
          mem_mb=50000,
          queue="medium"
      threads: 4
      script: 
          "02_permutation_other_rand_report.Rmd"

#-------------------------------------------------------------------------------

if RUN_MARK_RAND_OTHER_PERM:

  rule perm_mark_rand_other:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_paths = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mark_{fraction}", fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          resolution_df = RESOLUTION_OTHER,
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          iterations = config["values"]["03_sce_analysis"]["nr_cores"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          cons_level_use = "conserved_markers",
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str       
      conda:
          "../../envs/reclust_scores_perm_others.yml"
      output:
          perm_score_df = OUTPUT_DAT_OTHER + "/06_pmrk/perm_score_df_{dataset}"
      resources:
          mem_mb=200000,
          queue = "long"
      threads: config["values"]["03_sce_analysis"]["nr_cores"]
      script: 
          "02_scripts_other/05_permutation_other_rand.R"

  # export pvals: compare conserved markers to random background
  rule mark_vs_rand:
      input:
          orig_score_df_list_input = rules.reclustering_other_scores.output.score_df_list,
          perm_score_df_input = rules.perm_mark_rand_other.output.perm_score_df
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "conserved_markers",
          cons_level_use_alt = "seu_mark",
          comparison = "mark-vs-rand"
      output:
          pval_score_df_output = OUTPUT_DAT_OTHER + "/08_expp/mark-vs-rand_{dataset}"
      script:
          "02_scripts_other/08_export_pval.R"
          
  # visualise reclustering permutation
  rule permutation_report_mark_rand:
      input:
          seu_list = rules.reclustering_other.output,
          score_df_list = rules.reclustering_other_scores.output,
          perm_score_df = rules.perm_mark_rand_other.output
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "conserved_markers"
      resources:
          mem_mb=50000,
          queue = "medium"
      output:
          OUTPUT_REP + "/conserved_markers/perm_cons_rand_{dataset}.html"
      script: 
          "02_permutation_other_rand_report.Rmd"

if RUN_MMMS_RAND_OTHER_PERM:

  rule perm_mmms_rand_other:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_paths = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mmms_{fraction}", fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          resolution_df = RESOLUTION_OTHER,
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          iterations = config["values"]["03_sce_analysis"]["nr_cores"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          cons_level_use = "mmusall_markers",
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str       
      conda:
          "../../envs/reclust_scores_perm_others.yml"
      output:
          perm_score_df = OUTPUT_DAT_OTHER + "/07_pmms/perm_score_df_{dataset}"
      resources:
          mem_mb=200000,
          queue="long"
      threads: config["values"]["03_sce_analysis"]["nr_cores"]
      script: 
          "02_scripts_other/05_permutation_other_rand.R"

  # export pvals: compare conserved markers to random background
  rule mmms_vs_rand:
      input:
          orig_score_df_list_input = rules.reclustering_other_scores.output.score_df_list,
          perm_score_df_input = rules.perm_mmms_rand_other.output.perm_score_df
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "mmusall_markers",
          cons_level_use_alt = "seu_mmms",
          comparison = "mmms-vs-rand"
      output:
          pval_score_df_output = OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-rand_{dataset}"
      resources:
          mem_mb=5000,
          queue="medium"
      threads: 4
      script:
          "02_scripts_other/08_export_pval.R"
          
#-------------------------------------------------------------------------------

if RUN_GNST_SIGN_RAND_OTHER_PERM:
  
  # requires conda channel genomedk
  rule perm_genesets_signature_rand_other:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_sign = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_sign_{fraction}", fraction = fractions),
          ensembl_mark = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mark_{fraction}", fraction = fractions),
          ensembl_mmms = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mmms_{fraction}", fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          resolution_df = RESOLUTION_OTHER,
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          iterations = config["values"]["03_sce_analysis"]["iterations"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str     
      conda:
          "../../envs/reclust_scores_perm_others.yml"
      output:
          perm_score_df_mark = OUTPUT_DAT_OTHER + "/07_perg/perm_score_df_mark_{dataset}",
          perm_score_df_mmms = OUTPUT_DAT_OTHER + "/07_perg/perm_score_df_mmms_{dataset}",
          perm_score_df_mmms_mark = OUTPUT_DAT_OTHER + "/07_perg/perm_score_df_mmms_mark_{dataset}"
      resources:
          mem_mb=200000,
          queue="long"
      threads: config["values"]["03_sce_analysis"]["nr_cores"]
      script: 
          "02_scripts_other/07_permutation_other_signrand.R"

  # visualise reclustering permutation
  rule permutation_genesets_report:
      input:
          score_df_list = rules.reclustering_other_scores.output,
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          perm_score_df_mark = rules.perm_genesets_signature_rand_other.output.perm_score_df_mark,
          perm_score_df_mmms = rules.perm_genesets_signature_rand_other.output.perm_score_df_mmms,
          perm_score_df_mmms_mark = rules.perm_genesets_signature_rand_other.output.perm_score_df_mmms_mark,
          ensembl_sign = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_sign_{fraction}", fraction = fractions),
          ensembl_mark = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mark_{fraction}", fraction = fractions),
          ensembl_mmms = expand(OUTPUT_DAT_OWN + "/02_endf/ensembl_mmms_{fraction}", fraction = fractions)
      params:
          resolution_df = RESOLUTION_OTHER,
          datasets_other_hsc = datasets_other_hsc,
          datasets_other_str = datasets_other_str
      resources:
          mem_mb=5000,
          queue="medium"
      output:
          OUTPUT_REP + "/genesets/perm_genesets_{dataset}.html"
      script: 
          "02_permutation_other_genesets_report.Rmd"
        
  # export pvals: compare conserved markers to signature + random
  rule mark_vs_signrand:
      input:
          # list of score dfs from all conservation levels and resolutions
          orig_score_df_list_input = rules.reclustering_other_scores.output.score_df_list,
          # df of permuted scores (from comparing gene sets) 
          perm_score_df_input = rules.perm_genesets_signature_rand_other.output.perm_score_df_mark
      params:
          resolution_df = RESOLUTION_OTHER, # determine correct resolution
          cons_level_use = "conserved_markers", # determine correct cons level
          cons_level_use_alt = "seu_mark", # alternative name for cons level
          comparison = "mark-vs-signrand" # record comparison
      resources:
          mem_mb=5000,
          queue="medium"
      output:
          pval_score_df_output = OUTPUT_DAT_OTHER + "/08_expp/mark-vs-signrand_{dataset}"
      script:
          "02_scripts_other/08_export_pval.R"  

  # export pvals: compare mmusall_markers (BL6 markers) to signature + random
  rule mmms_vs_signrand:
      input:
          orig_score_df_list_input = rules.reclustering_other_scores.output.score_df_list,
          perm_score_df_input = rules.perm_genesets_signature_rand_other.output.perm_score_df_mmms
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "mmusall_markers",
          cons_level_use_alt = "seu_mmms",
          comparison = "mmms-vs-signrand"
      output:
          pval_score_df_output = OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-signrand_{dataset}"
      resources:
          mem_mb=5000,
          queue="medium"
      script:
          "02_scripts_other/08_export_pval.R" 

  """
  # export pvals: compare mmusall_markers (all BL6 markers) to conserved markers + random
  rule mmms_vs_markrand:
      input:
          orig_score_df_list_input = rules.reclustering_other_scores.output.score_df_list,
          perm_score_df_input = rules.permutation_genesets.output.perm_score_df_mmms_mark
      params:
          resolution_df = RESOLUTION_OTHER,
          cons_level_use = "mmusall_markers",
          cons_level_use_alt = "seu_mmms",
          comparison = "mmms-vs-markrand"
      output:
          pval_score_df_output = OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-markrand_{dataset}"
      script:
          "02_scripts_other/08_export_pval.R" 
  """        


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# correct all pvalues (to be used)

# only use and correct the pvals that WILL be used for the publication
# other tests will eventually be removed once decision is final
if RUN_PVAL_CORRECTION:
  
  # check this later
  rule pval_correction:
      input:
          own_sign_rand = expand(OUTPUT_DAT_OWN + "/08_expp/sign-vs-rand_{fraction}", fraction = fractions),
          own_mark_rand = expand(OUTPUT_DAT_OWN + "/08_expp/mark-vs-rand_{fraction}", fraction = fractions),
          own_mmms_rand = expand(OUTPUT_DAT_OWN + "/08_expp/mmms-vs-rand_{fraction}", fraction = fractions),
          own_mark_signrand = expand(OUTPUT_DAT_OWN + "/08_expp/mark-vs-signrand_{fraction}", fraction = fractions),
          own_mmms_signrand = expand(OUTPUT_DAT_OWN + "/08_expp/mmms-vs-signrand_{fraction}", fraction = fractions),
          other_sign_rand = expand(OUTPUT_DAT_OTHER + "/08_expp/sign-vs-rand_{dataset}", dataset = datasets_other),
          #other_mark_rand = expand(OUTPUT_DAT_OTHER + "/08_expp/mark-vs-rand_{dataset}", dataset = datasets_other),
          other_mmms_rand = expand(OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-rand_{dataset}", dataset = datasets_other),
          #other_mark_signrand = expand(OUTPUT_DAT_OTHER  + "/08_expp/mark-vs-signrand_{dataset}", dataset = datasets_other),
          other_mmms_signrand = expand(OUTPUT_DAT_OTHER + "/08_expp/mmms-vs-signrand_{dataset}", dataset = datasets_other)
      output:
          pval_corrected_df = OUTPUT_DAT_OTHER + "/09_crpv/all_corrected_pval"
      resources:
          mem_mb=5000,
          queue="medium"
      script:
          "02_scripts_other/09_correct_all_pvals.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
FINAL REPORTS

# visualise one chosen re-clustering of other datasets with the best combination
# of scores
# for random features, always the same resolution as the signature, 
# because each random set of genes (later for permutation) might have a dif-
# ferent optimal resolution. So choosing the same as the signature
# is the easiest option
"""


rule reclustering_other_report_final:
    input:
        seu_list = rules.reclustering_other.output,
        score_df_list = rules.reclustering_other_scores.output,
        resolution_df = RESOLUTION_OTHER,
        frame_df = FRAME_OTHER
    output:
        OUTPUT_REP + "/final/reclustering_other_final_report_{dataset}.html"
    resources:
        mem_mb=50000,
        queue="medium"
    threads: 4
    script: 
      "02_reclustering_other_report_final.Rmd"

