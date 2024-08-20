#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures"

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
COLORS = config["base"] + config["metadata_paths"]["colors"]

CELL_TYPES_EXCLUDE = config["values"]["03_sce_analysis"]["cell_types_exclude"]
RESOLUTION_OTHER = config["base"] + config["metadata_paths"]["resolution_other"]
print(CELL_TYPES_EXCLUDE)
print(RESOLUTION_OTHER)

RUN_PERM_TEST = config["run_permutation_test"]

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"]
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"]
ENSEMBL_ZEB = config["base"] + config["metadata_paths"]["ensembl_zeb"]
ENSEMBL_NMR = config["base"] + config["metadata_paths"]["ensembl_nmr"]

#RECLUSTER_OTHER = config["recluster_other"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
print(METADATA)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
datasets_other = ["ts_all_stromal", "ts_bone_marrow", "ts_hscs_progenitors", "li_all_stromal", "nmr_sorted_hspc", "zeb_all_hspc"]

#-------------------------------------------------------------------------------

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_sign/signature_list_" + f]

targets = targets + [OUTPUT_REP + "/signatures_summary.html"]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_sign_" + f]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mark_" + f]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_ndge_" + f]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mmms_" + f]

# reclustering our own datasets
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/03_rclo/sce_" + f]
  targets = targets + [OUTPUT_DAT + "/03_rclo/score_df_" + f]
  targets = targets + [OUTPUT_REP + "/reclustering_own/reclustering_own_report_" + f + ".html"]

# reclustering other datasets
for d in datasets_other:
  targets = targets + [OUTPUT_DAT + "/04_rcls/reclustered_" + d + "_list"]
  targets = targets + [OUTPUT_DAT + "/04_rcls/score_df_" + d + "_list"]
  targets = targets + [OUTPUT_REP + "/reclustering_other/reclustering_other_report_" + d + ".html"]
  targets = targets + [OUTPUT_REP + "/reclustering_other/reclustering_other_selected_report_" + d + ".html"]
 
  # testing reclustering scores
  #targets = targets + [OUTPUT_REP + "/reclustering_scores/test_reclustering_scores_" + d + ".html"]

  # permutation
  #if RUN_PERM_TEST:
    #targets = targets + [OUTPUT_DAT + "/05_perm/" + d + "_perm_score_df"]
    #targets = targets + [OUTPUT_REP + "/reclustering_permutation_report_" + r + ".html"]

 
#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""   
# EXTRACT SIGNATURES
# 
# Get conserved signatures = 
# genes that are conserved marker genes provided by Veronica Busa and that
# are also non-differentially expressed
"""

# export list of data on marker genes, conserved signatures, and nDGEs
rule export_signature:
    input:
        celltype_ndge_list = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/06_nres/PC_0.05_FC_1.5/shared_genes_{fraction}_celltypes",
        marker_cons = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_{fraction}.RData",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
    output:
        signature_list = OUTPUT_DAT + "/01_sign/signature_list_{fraction}"
    params:
        cts_exclude = CELL_TYPES_EXCLUDE
    script:
        "scripts/01_export_signatures.R"    

# visualise signatures across species and cell types
rule signature_summary:
    input: 
        signature_list_hsc = OUTPUT_DAT + "/01_sign/signature_list_hsc",
        signature_list_str = OUTPUT_DAT + "/01_sign/signature_list_str",
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/signatures_summary.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "signatures_summary.Rmd"

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
# - conserved marker genes
# - nDGEs
# - conserved signatures
# - all BL6 marker genes
#
# Additionally, permutation tests are performed using random gene sets of the 
# same number of genes for reclustering for the three other species
"""

#-------------------------------------------------------------------------------
# prepare ensembl conversion tables for each gene set to be tested
# ensembl datasets are downloaded in prepare_datasets_snakefile.py
rule prepare_ensembl:
    input:
        signature_list = rules.export_signature.output,
        sce_inp = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    output:
        ensembl_sign = OUTPUT_DAT + "/02_endf/ensembl_sign_{fraction}",
        ensembl_mark = OUTPUT_DAT + "/02_endf/ensembl_mark_{fraction}",
        ensembl_ndge = OUTPUT_DAT + "/02_endf/ensembl_ndge_{fraction}",
        ensembl_mmms = OUTPUT_DAT + "/02_endf/ensembl_mmms_{fraction}"
    script:
        "scripts/02_prepare_ensembl.R"
      
#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING OUR DATA
#
# Data will be reclustered exactly as the original dataset with 2000 HVGs 
# was clustered, based on batch-corrected PC coordinates.
# Here, SCE object are subsetted to each gene set, PCA is performed again using
# non-batch corrected values, then clustering is performed at identical
# resolution and k = number of neighbors using the original functions
"""
rule reclustering_own:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        signature_list = rules.export_signature.output
    params:
        k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
        resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
        cts_exclude = CELL_TYPES_EXCLUDE,
        nr_cores = config["values"]["03_sce_analysis"]["nr_cores"] 
    output:
        sce_output = OUTPUT_DAT + "/03_rclo/sce_{fraction}"
    script:
        "scripts/03_reclustering_own.R"

#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING OTHER DATA
#
# Datasets from other species are prepared for analysis using 
# prepare_datasets_snakefile.py
#
# Standard Seurat approach with standard options starting from raw counts
# but with aforementioned gene sets instead of HVGs.
"""
rule reclustering_other:
    input:
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/{dataset}",
        ensembl_sign = expand(rules.prepare_ensembl.output.ensembl_sign, fraction = fractions),
        ensembl_mark = expand(rules.prepare_ensembl.output.ensembl_mark, fraction = fractions),
        ensembl_mmms = expand(rules.prepare_ensembl.output.ensembl_mmms, fraction = fractions)
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        cut_off_counts = config["values"]["03_sce_analysis"]["reclustering_cutoff_counts"], 
        nr_cores = config["values"]["03_sce_analysis"]["nr_cores"] 
    output:
        seu_output = OUTPUT_DAT + "/04_rcls/reclustered_{dataset}_list"
    script: 
        "scripts/04_reclustering_other.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING SCORES
# 
# In order to test how well a certain gene set was able to re-cluster a 
# dataset, we test different established (and several new) metrics for
# comparing two clusterings, or for evaluating the quality of clustering.
# We use the re-clustered other datasets to test and evaluate each score.
# We choose to use following metrics for evaluating re-clustering:
#
# - adjusted rand index
# - fowlkes-mallows index
# - variation of information
# - mean cluster purity of reclustered labels
# - proportion of fields that are 0 (own score)
# - mean proportion of cells of a cell type per cluster (own score)
#
# We use many different scores because some of these scores increase or
# decrease with a higher clustering resolution/nr of clusters and we aim to 
# reduce this bias by using several different kinds of scores.
"""

rule test_reclustering_scores:
    input: 
        seu_list = OUTPUT_DAT + "/04_rcls/reclustered_{dataset}_list"
    output:
        OUTPUT_REP + "/reclustering_scores/test_reclustering_scores_{dataset}.html"
    conda:
        "../../envs/reclustering_scores.yml"
    script:
        "test_reclustering_scores.Rmd"
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
# REPORTS
"""

# get the reclustering scores for our own re-clustered datasets
rule reclustering_own_scores:
    input:
        sce_input = rules.reclustering_own.output,
    params:
        cts_exclude = CELL_TYPES_EXCLUDE,
        reclustering_functions = "../../source/sce_functions_reclustering.R",
    conda:
        "../../envs/reclustering_scores.yml"
    output:
        score_df = OUTPUT_DAT + "/03_rclo/score_df_{fraction}"
    script:
        "scripts/03_reclustering_own_scores.R"
        
# visualise re-clustering of own datasets, including scores
rule reclustering_own_report:
    input: 
        sce_input = rules.reclustering_own.output,
        score_df = rules.reclustering_own_scores.output
    output:
        OUTPUT_REP + "/reclustering_own/reclustering_own_report_{fraction}.html"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "reclustering_own_report.Rmd"

#-------------------------------------------------------------------------------  

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
        "scripts/04_reclustering_other_scores.R"
        
# visualise re-clustering of other datasets, including with scores
rule reclustering_other_report:
    input:
        seu_list = rules.reclustering_other.output,
        score_df_list = rules.reclustering_other_scores.output
    #params:
        #colors_path = COLORS,
        #colors = "../../source/colors.R"
    output:
        OUTPUT_REP + "/reclustering_other/reclustering_other_report_{dataset}.html"
    script: 
        "reclustering_other_report.Rmd"

# visualise one chosen re-clustering of other datasets with the best combination
# of scores
# for random features, always the same resolution as the conserved signature, 
# because each random set of genes (later for permutation) might have a dif-
# ferent optimal resolution. So choosing the same as the conserved signature
# is the easiest option
rule reclustering_other_report_selected:
    input:
        seu_list = rules.reclustering_other.output,
        score_df_list = rules.reclustering_other_scores.output,
        resolution_df = RESOLUTION_OTHER
    #params:
        #colors_path = COLORS,
        #colors = "../../source/colors.R"
    output:
        OUTPUT_REP + "/reclustering_other/reclustering_other_selected_report_{dataset}.html"
    script: 
        "reclustering_other_report_selected.Rmd"



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""
# PERMUTATION TEST

# Do permutation tests to see if using conserved signature gene sets is 
# significantly better than using the same number of random genes (as
# background).
# re-cluster seurat objects with random genes n = iteration times and compare
# the re-clustering scores of the original seu-object to the normal distribu-
# tion of the permuted scores.

"""
if RUN_PERM_TEST:

  # requires conda channel genomedk
  rule permutation_test:
      input:
          seu_preprocessed = config["base"] + config["metadata_paths"]["reclustering"] + "/prepared/{dataset}",
          ensembl_sign = expand(rules.prepare_ensembl.output.ensembl_sign, fraction = fractions)
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          iterations = 2,
          resolution_df = RESOLUTION_OTHER,
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          cut_off_counts = config["values"]["03_sce_analysis"]["reclustering_cutoff_counts"]       
      conda:
          "../../envs/reclustering_permutation.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/05_perm/{dataset}_perm_score_df"
      script: 
          "scripts/05_permutation_test.R"

"""

# visualise reclustering permutation
rule reclustering_permutation_report:
    input:
        seu_list = OUTPUT_DAT + "/04_rcls/reclustered_{reference}_list",
        score_df = rules.permutation_test.output.perm_score_df
    params:
        plotting = "../../source/plotting.R",
        resolution = resolution_list
    output:
        OUTPUT_REP + "/reclustering_permutation_report_{reference}.html"
    script: 
        "reclustering_permutation_report.Rmd"
"""

"""
# evaluate cluster silhuoette and purity of all re-clustering
rule reclustering_hum_eval:
    input:
        seu_list = OUTPUT_DAT + "/04_rcls/reclustered_{reference}_list"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    output:
        OUTPUT_REP + "/reclustering_hum_eval_{reference}.html"
    script: 
        "reclustering_hum_eval.Rmd"
"""

