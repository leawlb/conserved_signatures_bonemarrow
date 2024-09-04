#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]

OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/01_reclustering_own"

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
COLORS = config["base"] + config["metadata_paths"]["colors"]

CELL_TYPES_EXCLUDE = config["values"]["03_sce_analysis"]["cell_types_exclude"]

RUN_OWN_SIGN_PERM = config["run_own_sign_permutation"]
RUN_PERM_GENESETS = config["run_permutation_genesets"]

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"]
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"]
ENSEMBL_ZEB = config["base"] + config["metadata_paths"]["ensembl_zeb"]
ENSEMBL_NMR = config["base"] + config["metadata_paths"]["ensembl_nmr"]

#-------------------------------------------------------------------------------

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
print(METADATA)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")

#-------------------------------------------------------------------------------

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_gens/geneset_list_" + f]

  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_sign_" + f]
  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mark_" + f]
  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_ndge_" + f]
  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mmms_" + f]
  
  targets = targets + [OUTPUT_DAT + "/03_recl/sce_" + f]
  # targets = targets + [OUTPUT_DAT + "/04_rcls/score_df_" + f]
  # targets = targets + [OUTPUT_REP + "/reclustering_own_report_" + f + ".html"]

  if RUN_PERM_GENESETS:
    targets = targets + [OUTPUT_DAT + "/05_perg/perm_score_df_mark_" + f]
    targets = targets + [OUTPUT_DAT + "/05_perg/perm_score_df_mmms_" + f]
    targets = targets + [OUTPUT_DAT + "/05_perg/perm_score_df_mmms_mark_" + f]
    #targets = targets + [OUTPUT_REP + "/perm_genesets_" + f + ".html"]

  if RUN_OWN_SIGN_PERM:
    targets = targets + [OUTPUT_DAT + "/06_psig/perm_score_df_" + f]
    #targets = targets + [OUTPUT_REP + "/perm_conserved_signature_" + f + ".html"]
          
targets = targets + [OUTPUT_REP + "/genesets_summary.html"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""   
# EXTRACT GENE SETS
# 
# Get different gene sets:
# - conserved signature genes = conserved marker genes + nDGEs (SIGN)
# - conserved marker genes (MARK)
# - all BL6 marker genes (MMMS)
# - nDGEs
"""

# export list of gene sets
rule export_genesets:
    input:
        celltype_ndge_list = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/06_nres/PC_0.05_FC_1.5/shared_genes_{fraction}_celltypes",
        marker_cons = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_{fraction}.RData",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
    output:
        geneset_list = OUTPUT_DAT + "/01_gens/geneset_list_{fraction}"
    params:
        cts_exclude = CELL_TYPES_EXCLUDE
    script:
        "01_scripts_own/01_export_genesets.R"    

# visualise gene sets across species and cell types, especially SIGN
rule genesets_summary:
    input: 
        geneset_list_hsc = OUTPUT_DAT + "/01_gens/geneset_list_hsc",
        geneset_list_str = OUTPUT_DAT + "/01_gens/geneset_list_str",
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/genesets_summary.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "01_genesets_summary.Rmd"

# prepare ensembl conversion tables for each gene set to be tested
# ensembl datasets are downloaded in 00_prepare_datasets_snakefile.py
# includes species for other datasets in 02_signatures_reclustering_other
rule prepare_ensembl:
    input:
        geneset_list = rules.export_genesets.output,
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
        "01_scripts_own/02_prepare_ensembl.R"
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
       
"""
# RECLUSTERING AND PERMUTATION TESTS
# 
# Reclustering datasets to test the ability of the different signatures
# to capture cell identity.
# Using our own datasets with:
#
# - conserved signature genes
# - conserved marker genes
# - all BL6 marker genes 
# - nDGEs 

# Important: subclustering genes are removed in cell-type specific manner

# Additionally, permutation tests are performed using random genes 
# of a geneset + random genes as backgrounds
"""

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING
#
# Data will be reclustered as the original dataset with some adjustments.
# - SCE object are subsetted to each gene set
# - PCA is performed again using non-batch corrected values

# Then, clustering is performed at identical resolution and k = number of 
# neighbors as before.
"""
rule reclustering_own:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        geneset_list = rules.export_genesets.output
    params:
        k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
        resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
        nr_cores = 4, # 4 genesets total for reclustering
        cts_exclude = CELL_TYPES_EXCLUDE,
        functions_reclustering = "../../source/sce_functions_reclustering.R"
    output:
        sce_output = OUTPUT_DAT + "/03_recl/sce_{fraction}"
        # contains re-clustered vectors
        # cts_exclude already removed
        # HSC sub-sampled to 25,000 cells
    script:
        "01_scripts_own/03_reclustering_own.R"

# get the reclustering scores for our own re-clustered datasets
rule reclustering_own_scores:
    input:
        sce_input = rules.reclustering_own.output,
    params:
        cts_exclude = CELL_TYPES_EXCLUDE,
        functions_reclustering = "../../source/sce_functions_reclustering.R"
    conda:
        "../../envs/reclustering_scores.yml"
    output:
        score_df = OUTPUT_DAT + "/04_rcls/score_df_{fraction}"
    script:
        "01_scripts_own/04_reclustering_own_scores.R"
        
# visualise re-clustering of own datasets, including scores
rule reclustering_own_report:
    input: 
        sce_input = rules.reclustering_own.output,
        score_df = rules.reclustering_own_scores.output
    output:
        OUTPUT_REP + "/reclustering_own_report_{fraction}.html"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "01_reclustering_own_report.Rmd"

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
# or conserved markers genes + x random genes

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

# permutation: diff marker genes vs. conserved signature + random (same number)
if RUN_PERM_GENESETS:

  rule permutation_genesets:
      input:
          sce_input = rules.reclustering_own.output,
          # cts_exclude already excluded
          # HSCs subsampled to 25,000 cells
          geneset_list = rules.export_genesets.output,
      params:
          functions_reclustering = "../../source/sce_functions_reclustering.R",
          k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          #nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          #iterations = config["values"]["03_sce_analysis"]["iterations"],
          iterations = 40,
          nr_cores = 20,
          cut_off_counts = config["values"]["03_sce_analysis"]["reclustering_cutoff_counts"],
          cts_exclude = CELL_TYPES_EXCLUDE,
      conda:
          "../../envs/reclustering_own_scores.yml"
      output:
          perm_score_df_mark = OUTPUT_DAT + "/05_perg/perm_score_df_mark_{fraction}",
          perm_score_df_mmms = OUTPUT_DAT + "/05_perg/perm_score_df_mmms_{fraction}",
          perm_score_df_mmms_mark = OUTPUT_DAT + "/05_perg/perm_score_df_mmms_mark_{fraction}"
      script: 
          "01_scripts_own/05_permutation_genesets.R"

  # visualise reclustering permutation
  rule permutation_genesets_report:
      input:
          score_df = rules.reclustering_own_scores.output,
          geneset_list = rules.export_genesets.output,
          perm_score_df_mark = rules.permutation_genesets.output.perm_score_df_mark,
          perm_score_df_mmms = rules.permutation_genesets.output.perm_score_df_mmms,
          perm_score_df_mmms_mark = rules.permutation_genesets.output.perm_score_df_mmms_mark
      output:
          OUTPUT_REP + "/perm_genesets_{fraction}.html"
      script: 
          "01_permutation_own_genesets_report.Rmd"

     
#-------------------------------------------------------------------------------

# permutation: conserved signature vs. random (same number)
if RUN_OWN_SIGN_PERM:

  rule permutation_own_background_sign:
      input:
          sce_input = rules.reclustering_own.output,
          # cts_exclude already excluded
          # HSCs subsampled to 25,000 cells
          geneset_list = rules.export_genesets.output,
      params:
          k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          cut_off_counts = config["values"]["03_sce_analysis"]["reclustering_cutoff_counts"],
          #iterations = config["values"]["03_sce_analysis"]["iterations"],
          #nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          iterations = 40,
          nr_cores = 20,
          functions_reclustering = "../../source/sce_functions_reclustering.R",
          cts_exclude = CELL_TYPES_EXCLUDE,       
          cons_level_use = "conserved_signature"      
      conda:
          "../../envs/reclustering_own_scores.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/06_psig/perm_score_df_{fraction}"
      script: 
          "01_scripts_own/06_permutation_background.R"

  # visualise conserved signature permutation
  rule permutation_own_background_sign_report:
      input:
          score_df = rules.reclustering_own_scores.output,
          perm_score_df = rules.permutation_own_background_sign.output
      params:
          cons_level_use = "conserved_signature"
      output:
          OUTPUT_REP + "/perm_conserved_signature_{fraction}.html"
      script: 
          "01_permutation_own_background_report.Rmd"

