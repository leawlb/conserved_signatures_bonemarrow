#!/bin/python 

# import pandas as pd

# the purpose of this snakefile is to
# - prepare publically available scRNAseq datasets and for re-clustering
# - download ensembl ID conversion tables between MMUS and other species

# all of the generated data is saved in METADATA

# last run 2024-08-06

#-------------------------------------------------------------------------------

DIR_RECLUSTERING = config["base"] + config["metadata_paths"]["reclustering"]

RAN_LI_IPYNB = config["ran_li_ipynb"]

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"]
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"]
ENSEMBL_ZEB = config["base"] + config["metadata_paths"]["ensembl_zeb"]
ENSEMBL_NMR = config["base"] + config["metadata_paths"]["ensembl_nmr"]

#-------------------------------------------------------------------------------

targets = []

# ensembl conversion sheets
targets = targets + [ENSEMBL_MUS]
targets = targets + [ENSEMBL_HUM]
targets = targets + [ENSEMBL_ZEB]
targets = targets + [ENSEMBL_NMR]

# downloaded datasets (tabula sapiens)
targets = targets + [DIR_RECLUSTERING + "/raw/ts_bone_marrow"]
targets = targets + [DIR_RECLUSTERING + "/raw/ts_hscs_progenitors"]
targets = targets + [DIR_RECLUSTERING + "/raw/ts_all_stromal"]

# pre-processed datasets (tabula sapiens)
targets = targets + [DIR_RECLUSTERING + "/pre-processed/ts_bone_marrow"]
targets = targets + [DIR_RECLUSTERING + "/pre-processed/ts_hscs_progenitors"]
targets = targets + [DIR_RECLUSTERING + "/pre-processed/ts_all_stromal"]

# pre-process dataset (Li) only after ipynb script was executed
if RAN_LI_IPYNB:
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_only_counts"]
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_nocounts"]
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_all"]
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_stromal"]
  targets = targets + [DIR_RECLUSTERING + "/pre-processed/li_all_stromal"]

#-------------------------------------------------------------------------------

wildcard_constraints: 
    fraction="[a-z]+"

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------

# download ensembl conversion tables 
rule download_ensembl:
    output:
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    script:
        "prepare_datasets/download_ensembl.R"

#-------------------------------------------------------------------------------
"""
Human
"""

# download datasets
# the directories "/Li/download" and "/Li/data" must be generated manually
rule download_datasets:
    output:
        tabula_sapiens_bone_marrow = DIR_RECLUSTERING + "/raw/ts_bone_marrow",
        tabula_sapiens_hsc_progenitors = DIR_RECLUSTERING + "/raw/ts_hscs_progenitors",
        tabula_sapiens_stromal = DIR_RECLUSTERING + "/raw/ts_all_stromal"
    params:
        li_bone_marrow = DIR_RECLUSTERING + "/Li"
    script:
        "prepare_datasets/download_files.R"
        
# prepare datasets - Tabula Sapiens
rule prepare_ts_datasets:
    input:
        ts_bone_marrow_input = rules.download_datasets.output.tabula_sapiens_bone_marrow,
        ts_hsc_progenitors_input = rules.download_datasets.output.tabula_sapiens_hsc_progenitors,
        ts_stromal_input = rules.download_datasets.output.tabula_sapiens_stromal
    output:
        ts_bone_marrow_output = DIR_RECLUSTERING + "/pre-processed/ts_bone_marrow",
        ts_hsc_progenitors_output = DIR_RECLUSTERING + "/pre-processed/ts_hscs_progenitors",
        ts_stromal_output = DIR_RECLUSTERING + "/pre-processed/ts_all_stromal"
    script:
        "prepare_datasets/prepare_ts_datasets.R"

"""
# Manually proprocess Li et al. human stromal bone marrow dataset
# https://doi.org/10.7554/elife.81656
# according to https://github.com/Hongzhe2022/MSC_BM_scripts.

# After manual pre-processing, convert to Seurat and prepared for analysis.

# Only execute this rule when Load_GEO_data_analyze_it.ipynb has been executed.
"""
if RAN_LI_IPYNB:
  rule preprocess_li:
      input:
          h5ad_input = DIR_RECLUSTERING + "/Li/download/GSE190965_analyzed.h5ad",
          h5ad_raw_input = DIR_RECLUSTERING + "/Li/download/GSE190965_raw.h5ad"
      output:
          sce_raw_output_1 = DIR_RECLUSTERING + "/Li/sce_li_only_counts",
          sce_output_1 = DIR_RECLUSTERING + "/Li/sce_li_nocounts",
          sce_output_2 = DIR_RECLUSTERING + "/Li/sce_li_all",
          sce_output_3 = DIR_RECLUSTERING + "/Li/sce_li_stromal",
          seurat_output = DIR_RECLUSTERING + "/pre-processed/li_all_stromal"
      conda:
          "../../envs/zellkonverter_li.yml"
      script:
          "prepare_datasets/convert_lidata.R"
        
#-------------------------------------------------------------------------------
"""
NMR & Zebrafish
"""
"""
rule download_datasets_other:
    output:
    script:
        "prepare_datasets/download_ensembl.R"  

"""
