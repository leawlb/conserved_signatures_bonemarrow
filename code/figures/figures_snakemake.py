#!/bin/python 

#-------------------------------------------------------------------------------

# generate figures 

#-------------------------------------------------------------------------------

import pandas as pd

#-------------------------------------------------------------------------------


# snakemake --cluster "bsub -q medium -n1 -R rusage[mem=40GB]" -p -j1 -c1 -s figures_snakemake.py --latency-wait 20 --configfile ../config.yaml -n


#-------------------------------------------------------------------------------

# base path
OUTPUT_BASE = config["base"]
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
# Figure 1

targets = targets + [OUTPUT_REP + "/figure1.html"]


#-------------------------------------------------------------------------------

# local execution of non-demanding rules
localrules: all, figure1

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
      sce_input_list = expand(OUTPUT_BASE + "/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10", fraction = fractions)
  output:
      OUTPUT_REP + "/figure1.html"
  params:
      colors_path = COLORS,
      colors = "../source/colors.R"
  script:
      "figure1.Rmd"











