#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/05_cellrank"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/05_cellrank"

COLORS = config["base"] + config["metadata_paths"]["colors"]

print(OUTPUT_DAT)
#-------------------------------------------------------------------------------

targets = []
targets = targets + [OUTPUT_DAT + "/01_adat/adata_hsc_01.h5ad"]
targets = targets + [OUTPUT_DAT + "/02_scnp/adata_hsc_02.h5ad"]
targets = targets + [OUTPUT_DAT + "/03_cllr/adata_hsc_03.h5ad"]
targets = targets + [OUTPUT_REP + "/python_plots.pdf"]
targets = targets + [OUTPUT_DAT + "/04_psce/sce_hsc_pseudotime"]
targets = targets + [OUTPUT_REP + "/cellrank_report.html"]
targets = targets + [OUTPUT_DAT + "/05_bsce/sce_ery"]
targets = targets + [OUTPUT_DAT + "/05_bsce/sce_lym"]
targets = targets + [OUTPUT_DAT + "/05_bsce/sce_neu"]

#-------------------------------------------------------------------------------

localrules: all  
rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
# convert and export HSC SCE object into a scanpy adata.h5ad object
# remove unneccessary data before conversion

# only hsc fraction will be used for pseudotime calculation
# cellrank.yml also has ipykernel installed to allow use for .ipynb notebooks
rule export_to_adata:
    input: 
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10"
    output:
        adata_output = OUTPUT_DAT + "/01_adat/adata_hsc_01.h5ad"
    conda:
        "../../envs/zellkonverter_export.yml"
    resources:
        mem_mb=30000,
        queue="medium"
    threads: 4
    script:
        "scripts/01_export_to_adata.R"


# basic scanpy pre-processing steps to prepare for cellrank
# calculate neighbors, diffusion, pseudotime
rule scanpy_preprocessing:
    input: 
        adata_input = rules.export_to_adata.output.adata_output
    output:
        adata_output = OUTPUT_DAT + "/02_scnp/adata_hsc_02.h5ad"
    conda:
        "../../envs/cellrank.yml"
    resources:
        mem_mb=10000,
        queue="medium"
    threads: 10
    script:
        "scripts/02_scanpy_preprocessing.py"
  
# use cellrank to define terminal states and calculate cell fate probabilities  
# export all plots
rule cellrank:
    input: 
        adata_input = rules.scanpy_preprocessing.output.adata_output
    output:
        adata_output = OUTPUT_DAT + "/03_cllr/adata_hsc_03.h5ad",
        pdf_output = OUTPUT_REP + "/python_plots.pdf"
    params:
        n_states_total = 6,
        n_states_terminal = 3
    conda:
        "../../envs/cellrank.yml"
    resources:
        mem_mb=20000,
        queue="medium"
    threads: 10
    script:
        "scripts/03_cellrank.py"

# convert back to a SCE object with relevant metadata and coldata
# including pseudotime, fate probabilities
rule import_to_sce:
    input: 
        adata_input = rules.cellrank.output.adata_output
    output:
        sce_output = OUTPUT_DAT + "/04_psce/sce_hsc_pseudotime"
    conda:
        "../../envs/zellkonverter_import.yml"
    resources:
        mem_mb=30000,
        queue="medium"
    threads: 4
    script:
        "scripts/04_import_to_sce.R"

#-------------------------------------------------------------------------------
rule cellrank_report:
    input: 
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_pseudotime = rules.import_to_sce.output.sce_output
    output:
        OUTPUT_REP + "/cellrank_report.html"
    params:
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    resources:
        mem_mb=100000,
        queue="medium"
    threads: 4
    script:
        "cellrank_report.Rmd"
        
#-------------------------------------------------------------------------------

# separate SCE into three objects, one for each branch
# SCE is separated based on the cell fate probability for each lineage/branch
rule separate_by_branch:
    input: 
        sce_pseudotime = rules.import_to_sce.output.sce_output,
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10"
    output:
        sce_ery = OUTPUT_DAT + "/05_bsce/sce_ery",
        sce_lym = OUTPUT_DAT + "/05_bsce/sce_lym",
        sce_neu = OUTPUT_DAT + "/05_bsce/sce_neu"
    resources:
        mem_mb=100000,
        queue="medium"
    threads: 4
    script:
        "scripts/05_separate_by_branch.R"


