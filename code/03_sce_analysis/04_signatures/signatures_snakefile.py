#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures"

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
COLORS = config["base"] + config["metadata_paths"]["colors"]

GENE_LIST_SUBCLUSTERING = config["base"] + config["metadata_paths"]["gene_list_subclustering"]
CELL_TYPES_EXCLUDE = config["values"]["03_sce_analysis"]["cell_types_exclude"]
print(CELL_TYPES_EXCLUDE)

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"],
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"],
ENSEMBL_ZEB = config["base"] + config["metadata_paths"]["ensembl_zeb"],
ENSEMBL_NMR = config["base"] + config["metadata_paths"]["ensembl_nmr"],

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
print(METADATA)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
references_human = ["ts_all_stromal", "ts_bone_marrow", "ts_hsc_progenitors", "li_all_stromal"]

#-------------------------------------------------------------------------------

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_sign/signature_list_" + f]
  targets = targets + [OUTPUT_DAT + "/02_rcls/sce_" + f]
  targets = targets + [OUTPUT_REP + "/clustering_own_report_" + f + ".html"] 
# # 
# targets = targets + [ENSEMBL_MUS]
# targets = targets + [ENSEMBL_HUM]
# targets = targets + [ENSEMBL_ZEB]
# targets = targets + [ENSEMBL_NMR]
# 
# targets = targets + [OUTPUT_DAT + "/04_endf/ensembl_sign_" + f]
# targets = targets + [OUTPUT_DAT + "/04_endf/ensembl_mark_" + f]
# targets = targets + [OUTPUT_DAT + "/04_endf/ensembl_ndge_" + f]

# for r in references_human:
#   targets = targets + [OUTPUT_DAT + "/04_rcls/reclustered_" + r + "_list"]
#   targets = targets + [OUTPUT_REP + "/reclustering_hum_eval_" + r + ".html"]
#   targets = targets + [OUTPUT_REP + "/reclustering_hum_report_" + r + ".html"]
#   targets = targets + [OUTPUT_DAT + "/05_perm/" + r + "_score_df"]
#   targets = targets + [OUTPUT_REP + "/reclustering_permutation_report_" + r + ".html"]
# 
#targets = targets + [OUTPUT_REP + "/signatures_summary.html"]

#-------------------------------------------------------------------------------

wildcard_constraints: 
    fraction="[a-z]+"

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
"""   
# Get overlap with genes with conserved marker gene function (cons_markers)
# provided by Veronica Busa

"""

# export list of data on marker genes, conserved signatures, and nDGEs
rule export_signature:
    input:
        celltype_ndge_list = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/06_nres/PC_0.05_FC_1.5/shared_genes_{fraction}_celltypes",
        marker_cons = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_{fraction}.RData",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10"
    output:
        signature_list = OUTPUT_DAT + "/01_sign/signature_list_{fraction}"
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
# Reclustering datasets to test the ability of our conserved signatures
# to capture cell identity
# our own datasets and from other species (human, zebrafish, naked mole rat)
# using:
#
# - conserved marker genes
# - nDGEs
# - conserved signatures
#
# Additionally, perform permutation tests using random gene sets of the 
# same number of genes for reclustering


# download ensembl conversion tables 
# this is metadata
rule download_ensembl:
    output:
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    script:
        "scripts/03_download_ensembl.R"

# prepare ensembl conversion tables for each gene set to be tested
rule prepare_ensembl:
    input:
        signature_list = rules.export_signature.output,
        sce_inp = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    output:
        ensembl_sign = OUTPUT_DAT + "/04_endf/ensembl_sign_{fraction}",
        ensembl_mark = OUTPUT_DAT + "/04_endf/ensembl_mark_{fraction}",
        ensembl_ndge = OUTPUT_DAT + "/04_endf/ensembl_ndge_{fraction}"
    script:
        "scripts/04_prepare_ensembl.R"
        
"""

#-------------------------------------------------------------------------------
"""
# Re-clustering our own data

"""
# reclustering using original clustering pipeline (scater/igraph)
rule reclustering_own:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        signature_list = rules.export_signature.output
    params:
        k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
        resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
        cts_exclude = CELL_TYPES_EXCLUDE
    output:
        sce_output = OUTPUT_DAT + "/02_rcls/sce_{fraction}"
    script:
        "scripts/02_reclustering_own.R"


# visualise reclustered datasets
rule reclustering_own_report:
    input: 
        sce_input = rules.reclustering_own.output
    output:
        OUTPUT_REP + "/clustering_own_report_{fraction}.html"
    params:
        colors_path = COLORS,
        #functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "reclustering_own_report.Rmd"



#-------------------------------------------------------------------------------
# reclustering test dataset 1 ("All Stromal" from tabula sapiens)
# standard Seurat approach with standard options from raw counts
# but without HVGs
"""
rule reclustering_human:
    input:
        seu_input = config["base"] + config["metadata_paths"]["human_test_datasets"] + "/{reference}",
        ensembl_emfs = expand(rules.prepare_ensembl.output.ensembl_emfs, fraction = fractions),
        ensembl_mark = expand(rules.prepare_ensembl.output.ensembl_mark, fraction = fractions)
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        cut_off_counts = 10,
        nr_cores = 20
    output:
        seu_output = OUTPUT_DAT + "/04_rcls/reclustered_{reference}_list"
    script: 
        "scripts/04_reclustering_hum.R"

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

# visualise reclustered datasets using different gene sets
rule reclustering_hum_report:
    input:
        seu_list = OUTPUT_DAT + "/04_rcls/reclustered_{reference}_list"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    output:
        OUTPUT_REP + "/reclustering_hum_report_{reference}.html"
    script: 
        "reclustering_hum_report.Rmd"

resolution_list = [
  {"li_all_stromal": 0.4},
  {"ts_all_stromal": 0.7},
  {"ts_hsc_progenitors": 0.35},
  {"ts_bone_marrow": 0.35}
]

print(resolution_list)

#-------------------------------------------------------------------------------
# perform permutation tests for random gene sets instead of 
# EMF sets to obtain a pval

# requires conda channel genomedk
rule permutation_test:
    input:
        seu_preprocessed = config["base"] + config["metadata_paths"]["human_test_datasets"] + "/{reference}",
        ensembl_emfs = expand(rules.prepare_ensembl.output.ensembl_emfs, fraction = fractions)
    params:
        reclustering_functions = "../../source/sce_functions_reclustering.R",
        iterations = 2000,
        nr_cores = 20,
        cut_off_counts = 10,
        resolution = resolution_list
    conda:
        "../../envs/reclustering_permutation.yml"
    output:
        perm_score_df = OUTPUT_DAT + "/05_perm/{reference}_score_df"
    script: 
        "scripts/05_permutation_test.R"

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
