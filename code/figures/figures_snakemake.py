#!/bin/python 

#-------------------------------------------------------------------------------

# generate figures 

#-------------------------------------------------------------------------------

import pandas as pd

#-------------------------------------------------------------------------------

# use this command in this repository with snkmk_isbm
"""
snakemake --cluster "bsub -q medium -n1 -R rusage[mem=80GB]" -p -j2 -c1 -s figures_snakemake.py --latency-wait 20 --configfile ../config.yaml -n
"""

#-------------------------------------------------------------------------------

# base path
BASE = config["base"]
OUTPUT_REP = config["base"] + "/data/manuscript1/figures"
print(OUTPUT_REP)

#-------------------------------------------------------------------------------

# load metadata table to obtain species, individuals, fractions to be used
METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])

# function to extract the required variables
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

# get required variables
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

#-------------------------------------------------------------------------------

# get global paths
COLORS = config["base"] + config["metadata_paths"]["colors"]

#-------------------------------------------------------------------------------

# define targets for rule all
targets = []

#-------------------------------------------------------------------------------

targets = targets + [OUTPUT_REP + "/figure1.html"]
targets = targets + [OUTPUT_REP + "/figure2.html"]


#-------------------------------------------------------------------------------

# local execution of non-demanding rules
localrules: all

# define rules
rule all: # must contain all possible output paths from all rules below and must always be the first rule
    input:
        targets
        

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# .Rmd files are generated so that they can both be used manually and run
# from this snakefile

"""
Figure 1
"""

rule figure1:
  input:
      # fully annotated SCE objects, one per fraction (02_/03_/)
      sce_input_list = expand(BASE + "/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10", fraction = fractions),
      # list of genes (in sequence) for dotplots
      gene_list_dtplt = config["base"] + config["metadata_paths"]["gene_list_dotplot"]
  output:
      OUTPUT_REP + "/figure1.html"
  params:
      colors_path = COLORS,
      colors = "../source/colors.R"
  script:
      "figure1.Rmd"



"""
Figure 2
"""
# list of genes (in sequence) for dotplots
GENE_LIST_DOTPLOT = config["base"] + config["metadata_paths"]["gene_list_dotplot"]

rule figure2:
  input:
      # fully annotated SCE object (HSPCs)
      sce_hsc = BASE + "/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
      # fully annotated SCE object (HSPCs) with diffusion pseudotime
      sce_pt = BASE + "/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/05_cellrank/04_psce/sce_hsc_pseudotime",
      # list of gene sets for HSPCs (conserved signature, conserved markers, etc.)
      geneset_list_hsc = BASE + "/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own/01_gens/geneset_list_hsc"
  output:
      OUTPUT_REP + "/figure2.html"
  params:
      colors_path = COLORS,
      colors = "../source/colors.R"
  script:
      "figure2.Rmd"








