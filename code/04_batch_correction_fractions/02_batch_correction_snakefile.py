#!/bin/python 

"""
Info in batch correction methods from:
Luecken et al. "Benchmarking atlas-level data integration in single-cell 
genomics", Nat Met 2022
Tran, Ang, Chevrier, Zhang et al. "A benchmark of batch effect correction 
methods for single-cell RNA sequencing data", Genome Biology 2020
"""
#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import sys 

# paths from config
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["raw"])
BATCH_USE = config["values"]["batch_correction"]["batch_use_fractions"] # which SCE col to use as batch

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")

# construct paths for all possible outputs/targets, required for rule all
targets = []
for f in fractions:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/08_rnrm_fractions/sce_" + f + "_" + BATCH_USE + "-08"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/09_mnnc_fractions/sce_" + f + "_" + BATCH_USE +"-09"]
  targets = targets + [OUTPUT_BASE_PATH + "/reports/04_batch_correction_fractions/" + f + "/batch_correction_fraction_report_" + f + "_" + BATCH_USE + ".html"]
  if BATCH_USE != "Object_ID":
    targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/09_srt3_fractions/sce_" + f + "_" + BATCH_USE +"-09"]
    
targets = targets + [OUTPUT_BASE_PATH + "/reports/04_batch_correction_fractions/batch_correction_fraction_summary_" + BATCH_USE + ".html"]
  
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
 
"""
Quick and basic renormalization.
Scaling requires identical composition in all batches and worsens
bioconservation so scaling will not be used by default
"""

rule renormalize:
    input:
        sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_{fraction}-07"
    output:
        sce_08 = OUTPUT_BASE_PATH + "/sce_objects/08_rnrm_fractions/sce_{fraction}_" + BATCH_USE + "-08"
    params:
        batch_use = BATCH_USE,
        rescale = config["rescale_for_batch_correction"]
    script:
        "scripts/08_renormalize.R"

#-------------------------------------------------------------------------------
    
"""
batch correct using MNNcorrect:

MNNcorrect is good at recovering DEGs from batch corrected data and 
bioconservation, but slow and does not perform well on batch correction 
Requires shared cell types between batches but no labels
(fastMNN) seems to balance batch effect removal and bioconservation
HVG selection improves performance but restricts analysis
Seurat cannot integrate this many samples (Object_ID)
"""
rule run_mnncorrect:
    input:
        sce_08 = rules.renormalize.output 
    output:
        sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnnc_fractions/sce_{fraction}_" + BATCH_USE + "-09"
    params:
        batch_use = BATCH_USE,
        rescale = config["rescale_for_batch_correction"], # condition
        hvgs_for_batch_correction = config["hvgs_for_batch_correction"],# condition
        mnn_fast = config["mnn_use_fast"], # condition
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], # number
        sce_functions = "../source/sce_functions.R" # this is the working dir
    script:
        "scripts/09_mnncorrect.R"
        
"""
batch correct using Seurat3

Seurat3 is good at batch correction and among best for multiple batch integration
but not great at recovering DEGs from batch corrected data
unbalanced towards stronger batch effect removal, but successful at removing species batch effects
Requires shared cell types between batches but no labels, scaling little effect
It is easiest to use HVGs and Renormalize by seurat means.
"""
if BATCH_USE != "Object_ID":
  rule run_seurat3:
      input:
          sce_07 = rules.merge_datasets_species.output 
      output:
          sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_srt3_ractions/sce_{fraction}_" + BATCH_USE + SEURAT_RED + "-09"
      params:
          batch_use = BATCH_USE,
          nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], # number
          nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
          seurat_reduction = SEURAT_RED, # which reduction method to use
          sce_functions = "../source/sce_functions.R" # we are in the working dir
      script:
          "scripts/09_seurat.R"


#-------------------------------------------------------------------------------

"""
Make batch correction reports 
"""

if BATCH_USE == "Object_ID":
  sce_09_input = rules.run_mnncorrect.output,
if BATCH_USE != "Object_ID": # seurat cannot deal with Object_ID as batch
  sce_09_input = [rules.run_mnncorrect.output, rules.run_seurat3.output]

print(sce_09_input)

# takes >4h for hscs
rule make_reports:
    input:
        sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_{fraction}-07",
        sce_08 = rules.renormalize.output,
        sce_09_input = sce_09_input
    output:
        OUTPUT_BASE_PATH + "/reports/04_batch_correction_fractions/{fraction}/batch_correction_fraction_report_{fraction}_" + BATCH_USE + ".html"
    params:
        batch_use = BATCH_USE,
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        color_tables = TABLES_PATH
    script:
        "batch_correction_fraction_reports.Rmd" 
   
print(expand(rules.run_mnncorrect.output, fraction = fractions))
rule make_summary_report:
    input:
        sce_09 = expand(rules.run_mnncorrect.output, fraction = fractions)
    params:
        fractions = fractions
    output:
        OUTPUT_BASE_PATH + "/reports/04_batch_correction_fractions/batch_correction_fraction_summary_" + BATCH_USE + ".html"
    script:
        "batch_correction_summary.Rmd" 
