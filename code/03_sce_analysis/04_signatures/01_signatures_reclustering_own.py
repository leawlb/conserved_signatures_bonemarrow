#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures"

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
COLORS = config["base"] + config["metadata_paths"]["colors"]

CELL_TYPES_EXCLUDE = config["values"]["03_sce_analysis"]["cell_types_exclude"]
print(CELL_TYPES_EXCLUDE)

#RUN_PERM_GENESETS = config["run_permutation_genesets"]
#RUN_PERM_BACKGROUND_SIGN = config["run_permutation_background_sign"]
#RUN_PERM_BACKGROUND_MARK = config["run_permutation_background_mark"]
RUN_OWN_BACKGROUND_SIGN = config["run_own_permutation_background_sign"]

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
datasets_other_hsc = ["ts_hscs_progenitors", "ts_bone_marrow", "mus_weinreb_hspc", "mus_tm_bonemarrow", "nmr_sorted_hspc", "zeb_all_hspc"]
datasets_other_str = ["ts_all_stromal", "li_all_stromal", "mus_tik_stromal", "mus_bar_stromal"]
datasets_other = datasets_other_hsc + datasets_other_str
print(datasets_other)

#-------------------------------------------------------------------------------

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_sign/signature_list_" + f]

targets = targets + [OUTPUT_REP + "/signatures_summary.html"]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_sign_" + f]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mark_" + f]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_ndge_" + f]
targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mmms_" + f]

# reclustering own datasets
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/03_rclo/sce_" + f]
  targets = targets + [OUTPUT_DAT + "/03_rclo/score_df_" + f]
  targets = targets + [OUTPUT_REP + "/reclustering_own/reclustering_own_report_" + f + ".html"]
  
  # permutation of own datasets
  if RUN_OWN_BACKGROUND_SIGN:
     targets = targets + [OUTPUT_DAT + "/08_test/perm_score_df_" + f]
     

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
# genes that are conserved marker genes + non-differentially expressed
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
# Using our own datasets with:
#
# - conserved signatures
# - conserved marker genes
# - all BL6 marker genes 
# - nDGEs (own datasets only)
# - random genes (other datasets only)
# - all conserved markers - random genes (other datasets only)
# - all BL6 markers - random genes (other datasets only)
#
# Additionally, permutation tests are performed.
"""

#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING OUR DATA
#
# Data will be reclustered as the original dataset (2000 HVGs, corrected PC).
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
# REPORT RECLUSTERING OWN:
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

# Re-cluster own datasets with random genes n = iteration times and compare
# the re-clustering scores of the original seurat object to the normal 
# distribution of the permuted scores.

# For permutation, the same resolution as the original gene set is chosen, 
# because each random set of genes might have a different optimal resolution. 

"""
#-------------------------------------------------------------------------------

if RUN_OWN_BACKGROUND_SIGN:

  # requires conda channel genomedk
  rule permutation_own_background_sign:
      input:
          sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
          ensembl = rules.prepare_ensembl.output.ensembl_sign
      params:
          k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          cts_exclude = CELL_TYPES_EXCLUDE,       
          reclustering_functions = "../../source/sce_functions_reclustering.R",
          iterations = 10,
          resolution_df = RESOLUTION_OTHER,
          #nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          nr_cores = 5,
          cut_off_counts = config["values"]["03_sce_analysis"]["reclustering_cutoff_counts"],
          cons_level_use = "conserved_signature"      
      conda:
          "../../envs/reclustering_scores_perm_bioconductor.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/08_test/perm_score_df_{fraction}"
      script: 
          "scripts/test_own_permutation.R"



