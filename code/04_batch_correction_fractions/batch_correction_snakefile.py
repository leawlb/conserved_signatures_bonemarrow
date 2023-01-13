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
SEURAT_RED = config["values"]["batch_correction"]["seurat_reduction"] # which Seurat method to use 

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

# construct paths for all possible outputs/targets, required for rule all
targets =[OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_" + f + "-07" for f in fractions]
  
if config["run_mnncorrect"]: # to avoid starting bc before the merged objects are generated
  for f in fractions:
    targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/08_rnrm_fractions/sce_" + f + "_" + BATCH_USE + "-08"]
    targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/09_mnncorrect_fractions/sce_" + f + "_" + BATCH_USE +"-09"]
if config["run_seurat3"]:
    targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/09_seurat3_fractions/sce_" + f + "_" + BATCH_USE + SEURAT_RED +"-09"]
    targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction_fractions/" + f + "/batch_correction_fraction_report_" + f + "_" + BATCH_USE + SEURAT_RED + ".html"]
  
if config["run_batch_correction_summary"]:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction_fractions/batch_correction_fraction_summary_" + SEURAT_RED + ".html"]
  
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
        
"""
Merge all datasets of one fraction into two fraction dataset
Input data is located in 02_preprocessing/06_sglr/
    
Merge all datasets into one teo fraction datasets 
"""
rule merge_datasets_fractions:
    input: 
        sce_06 = expand(OUTPUT_BASE_PATH + "/sce_objects/06_sglr/{s}/", s = species)
    output:
        sce_07_hsc = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_hsc-07",
        sce_07_str = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_str-07"
    params:
        individuals = individuals,
        samples_to_remove = config["samples_to_remove"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
        species = species,
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    script:
        "scripts/07_merge_datasets.R" 


"""
Quick and basic renormalization.
HVG selection improves performance but restricts analysis
Scaling requires identical composition in all batches and worsens
bioconservation so scaling will not be used by default
"""
if config["run_mnncorrect"]:
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
"""
if config["run_mnncorrect"]:
  print("run_mnncorrect")
  rule run_mnncorrect:
      input:
          sce_08 = rules.renormalize.output 
      output:
          sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnncorrect_fractions/sce_{fraction}_" + BATCH_USE + "-09"
      params:
          batch_use = BATCH_USE,
          rescale = config["rescale_for_batch_correction"], # rule
          hvgs_for_batch_correction = config["hvgs_for_batch_correction"],# rule
          mnn_fast = config["mnn_use_fast"], # rule
          nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
          nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], # number
          sce_functions = "../source/sce_functions.R" # this is the working dir
      script:
          "scripts/09_mnncorrect.R"

"""
batch correct using Seurat3

Seurat3 CCA is good at batch correction and among best for multiple batch integration
but not great at recovering DEGs from batch corrected data
unbalanced towards stronger batch effect removal, but successful at removing species batch effects
Requires shared cell types between batches but no labels, scaling little effect
It is easiest to use HVGs and Renormalize by seurat means.

RPCA is better suited for imbalanced datasets and runs faster.
"""
if config["run_seurat3"]:
    rule run_seurat3:
        input:
            sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_{fraction}-07"
        output:
            sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3_fractions/sce_{fraction}_" + BATCH_USE + SEURAT_RED + "-09"
        params:
            batch_use = BATCH_USE,
            nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], # number
            nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
            seurat_reduction = SEURAT_RED, # which reduction method to use
            sce_functions = "../source/sce_functions.R" # this is the working dir
        script:
            "scripts/09_seurat.R"

#-------------------------------------------------------------------------------

"""
Make batch correction reports for each species and each method
due to multiple runs, different reports are stored in 04_correction_COLLECTION
"""
if config["run_seurat3"]:
  rule make_reports:
      input:
          sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_fractions/sce_{fraction}-07",
          sce_08 = rules.renormalize.output,
          sce_09mnn = rules.run_mnncorrect.output,
          sce_09seurat = rules.run_seurat3.output
      output:
          OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction_fractions/{fraction}/batch_correction_fraction_report_{fraction}_" + BATCH_USE + SEURAT_RED + ".html"
      params:
          nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
          species = species,
          color_tables = TABLES_PATH
      script:
          "batch_correction_fraction_reports.Rmd" 
        
if config["run_seurat3"]:
  rule make_summary_report:
      input:
          sce_09_path = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3"
      output:
          OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction_fractions/batch_correction_fraction_summary_" + SEURAT_RED + ".html"
      script:
          "batch_correction_summary.Rmd" 
