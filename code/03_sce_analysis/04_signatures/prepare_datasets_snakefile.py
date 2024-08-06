#!/bin/python 

# import pandas as pd

# the purpose of this snakefile is to
# - prepare publically available scRNAseq datasets and for re-clustering
# - download ensembl ID conversion tables between MMUS and other species

# all of the generated data is saved in METADATA

# last run 2024-08-06

#-------------------------------------------------------------------------------

DIR_RECLUSTERING = config["base"] + config["metadata_paths"]["reclustering"]

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
targets = targets + [DIR_RECLUSTERING + "/TEST/raw/ts_bone_marrow"]
targets = targets + [DIR_RECLUSTERING + "/TEST/raw/ts_hscs_progenitors"]
targets = targets + [DIR_RECLUSTERING + "/TEST/raw/ts_all_stromal"]

# pre-processed datasets (tabula sapiens)
targets = targets + [DIR_RECLUSTERING + "/TEST/pre-processed/ts_bone_marrow"]
targets = targets + [DIR_RECLUSTERING + "/TEST/pre-processed/ts_hscs_progenitors"]
targets = targets + [DIR_RECLUSTERING + "/TEST/pre-processed/ts_all_stromal"]

# prepare dataset (Li)
targets = targets + [DIR_RECLUSTERING + "/TEST/Li/download/GSE190965_analyzed.h5ad"]

# pre-process dataset (Li)
targets = targets + [DIR_RECLUSTERING + "/TEST/Li/sce_li_only_counts"]
targets = targets + [DIR_RECLUSTERING + "/TEST/Li/sce_li_nocounts"]
targets = targets + [DIR_RECLUSTERING + "/TEST/Li/sce_li_all"]
targets = targets + [DIR_RECLUSTERING + "/TEST/Li/sce_li_stromal"]
targets = targets + [DIR_RECLUSTERING + "/TEST/pre-processed/li_all_stromal"]



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
        tabula_sapiens_bone_marrow = DIR_RECLUSTERING + "/TEST/raw/ts_bone_marrow",
        tabula_sapiens_hsc_progenitors = DIR_RECLUSTERING + "/TEST/raw/ts_hscs_progenitors",
        tabula_sapiens_stromal = DIR_RECLUSTERING + "/TEST/raw/ts_all_stromal"
    params:
        li_bone_marrow = DIR_RECLUSTERING + "/TEST/Li"
    script:
        "prepare_datasets/download_files.R"
        
# prepare datasets - Tabula Sapiens
rule prepare_ts_datasets:
    input:
        ts_bone_marrow_input = rules.download_datasets.output.tabula_sapiens_bone_marrow,
        ts_hsc_progenitors_input = rules.download_datasets.output.tabula_sapiens_hsc_progenitors,
        ts_stromal_input = rules.download_datasets.output.tabula_sapiens_stromal
    output:
        ts_bone_marrow_output = DIR_RECLUSTERING + "/TEST/pre-processed/ts_bone_marrow",
        ts_hsc_progenitors_output = DIR_RECLUSTERING + "/TEST/pre-processed/ts_hscs_progenitors",
        ts_stromal_output = DIR_RECLUSTERING + "/TEST/pre-processed/ts_all_stromal"
    script:
        "prepare_datasets/prepare_ts_datasets.R"

"""
# proprocess Li et al. human stromal bone marrow dataset
# https://doi.org/10.7554/elife.81656
# according to https://github.com/Hongzhe2022/MSC_BM_scripts.
# 
# This requires flexible channel priority & reordering channels in the yml file:
# - bioconda
# - pkgs/main
# - conda-forge
# This is why not the original Load_GEO_data_analyze_it.yml in /Li/downloads but 
# Load_GEO_data_analyze_it_export.yml in /envs is used.
"""
rule prepare_li:
    output:
        # this is just the main output, but other intermediate files are also generated in the download folder
        DIR_RECLUSTERING + "/TEST/Li/download/GSE190965_analyzed.h5ad"
    conda:
        "../../envs/Load_GEO_data_analyze_it_export.yml"
    params:
        path_to_jupyter_notebook = DIR_RECLUSTERING + "/Li/download/Load_GEO_data_analyze_it.ipynb"
    shell:
        "jupyter nbconvert --execute {params.path_to_jupyter_notebook} "

# the pre-processed dataset is then used as input for conversion to Seurat
# and prepared for downstream analysis
rule preprocess_li:
    input:
        h5ad_input = DIR_RECLUSTERING + "/Li/download/GSE190965_analyzed.h5ad",
        h5ad_raw_input = DIR_RECLUSTERING + "/Li/download/GSE190965_raw.h5ad"
    output:
        sce_raw_output_1 = DIR_RECLUSTERING + "/TEST/Li/sce_li_only_counts",
        sce_output_1 = DIR_RECLUSTERING + "/TEST/Li/sce_li_nocounts",
        sce_output_2 = DIR_RECLUSTERING + "/TEST/Li/sce_li_all",
        sce_output_3 = DIR_RECLUSTERING + "/TEST/Li/sce_li_stromal",
        seurat_output = DIR_RECLUSTERING + "/TEST/pre-processed/li_all_stromal"
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
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    script:
        "prepare_datasets/download_ensembl.R"  

"""
