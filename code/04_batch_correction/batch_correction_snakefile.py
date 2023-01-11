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
BATCH_USE = config["values"]["batch_correction"]["batch_use"] # which SCE col to use as batch

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

species_build = species
if config["merge_all"]:
    species_build = species_build + ["all"]
print(species_build)

# construct paths for all possible outputs/targets, required for rule all
targets =[OUTPUT_BASE_PATH + "/sce_objects/07_mrge/sce_" + s + "-07" for s in species_build]
  
for s in species_build:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/08_rnrm/sce_" + s + "_" + BATCH_USE + "-08"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/09_mnncorrect/sce_" + s + "_" + BATCH_USE +"-09"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_" + s + "_" + BATCH_USE +"-09"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction/" + s + "/batch_correction_species_report_" + s + "_" + BATCH_USE + ".html"]
  
if config["run_batch_correction_summary"]:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction/batch_correction_species_summary.html"]
  
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
        
"""
Merge all datasets of one species into one species dataset
Input data is located in 02_preprocessing/04_norm/
"""

rule merge_datasets_species:
    input: 
        sce_06 = OUTPUT_BASE_PATH + "/sce_objects/06_sglr/{species}/"
    output:
        sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge/sce_{species}-07"
    params:
        individuals = individuals,
        samples_to_remove = config["samples_to_remove"],
        sce_functions = "../source/sce_functions.R", # this is the working dir
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"]
    wildcard_constraints:
        species = "[a-z]+"
    script:
        "scripts/07_merge_datasets.R"
      
"""       
Merge all datasets into one big dataset using the same script
"""
if config["merge_all"]:
  print("Merge all = True")
  rule merge_datasets_all:
      input: 
          sce_06 = expand(OUTPUT_BASE_PATH + "/sce_objects/06_sglr/{s}/", s = species)
      output:
          sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge/sce_all-07"
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

rule renormalize:
    input:
        sce_07 = rules.merge_datasets_species.output
    output:
        sce_08 = OUTPUT_BASE_PATH + "/sce_objects/08_rnrm/sce_{species}_" + BATCH_USE + "-08"
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
          sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_mnncorrect/sce_{species}_" + BATCH_USE + "-09"
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

Seurat3 is good at batch correction and among best for multiple batch integration
but not great at recovering DEGs from batch corrected data
unbalanced towards stronger batch effect removal, but successful at removing species batch effects
Requires shared cell types between batches but no labels, scaling little effect
It is easiest to use HVGs and Renormalize by seurat means.
"""
if config["run_seurat3"]:
    rule run_seurat3:
        input:
            sce_07 = rules.merge_datasets_species.output 
        output:
            sce_09 = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3/sce_{species}_" + BATCH_USE + "-09"
        params:
            batch_use = BATCH_USE,
            nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], # number
            nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
            sce_functions = "../source/sce_functions.R" # this is the working dir
        script:
            "scripts/09_seurat.R"

#-------------------------------------------------------------------------------

"""
Make batch correction reports for each species and each method
due to multiple runs, different reports are stored in 04_correction_COLLECTION
"""
rule make_reports:
    input:
        sce_07 = rules.merge_datasets_species.output,
        sce_08 = rules.renormalize.output,
        sce_09mnn = rules.run_mnncorrect.output,
        sce_09seurat = rules.run_seurat3.output
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction/{species}/batch_correction_species_report_{species}_" + BATCH_USE + ".html"
    params:
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        species_build = species_build,
        color_tables = TABLES_PATH
    script:
        "batch_correction_species_reports.Rmd" 
        
if config["run_batch_correction_summary"]:
  rule make_summary_report:
      input:
          sce_09_path = OUTPUT_BASE_PATH + "/sce_objects/09_seurat3"
      output:
          OUTPUT_BASE_PATH + "/sce_objects/reports/04_batch_correction/batch_correction_species_summary.html"
      script:
          "batch_correction_summary.Rmd" 
