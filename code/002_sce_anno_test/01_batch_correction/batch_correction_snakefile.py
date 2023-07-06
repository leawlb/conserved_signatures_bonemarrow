#!/bin/python 

"""
Info on batch correction methods from:
Luecken et al. "Benchmarking atlas-level data integration in single-cell 
genomics", Nat Met 2022
Tran, Ang, Chevrier, Zhang et al. "A benchmark of batch effect correction 
methods for single-cell RNA sequencing data", Genome Biology 2020
"""
#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
OUTPUT_BASE = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["table"])
BATCH_USE = config["values"]["batch_correction"]["batch_use_fractions"] # which colData to use as batch
# seurat kept for QC sake but not used for downstream analysis
SEURAT_RED = config["values"]["batch_correction"]["seurat_reduction"] # which seurat method to use

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno/"
OUTPUT_REP = OUTPUT_BASE + "/reports/02_sce_anno/01_batch_correction/"

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "01_mnnc/sce_" + f + "_" + BATCH_USE + "-01"]
  targets = targets + [OUTPUT_REP + "batch_correction_report_" + f + "_" + BATCH_USE + ".html"]
  if BATCH_USE != "Object_ID":
    targets = targets + [OUTPUT_DAT + "01_srt3/sce_" + f + "_" + BATCH_USE + "_" + SEURAT_RED +"-01"]
    
#targets = targets + [OUTPUT_REP + "batch_correction_summary_" + BATCH_USE + ".html"]
  
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
 
#-------------------------------------------------------------------------------
    
"""
batch correct using MNNcorrect:

MNNcorrect is good at recovering DEGs from batch corrected data and 
bioconservation, but slow and does not perform well on batch correction 
Requires shared cell types between batches but no labels
(fastMNN) seems to balance batch effect removal and bioconservation
HVG selection improves performance but restricts analysis
Seurat cannot integrate this many samples (Object_ID).

MNN performs its own normalization step.
"""
rule run_mnncorrect:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/fractions/sce_{fraction}-08"
    output:
        sce_output = OUTPUT_DAT + "01_mnnc/sce_{fraction}_" + BATCH_USE + "-01"
    params:
        batch_use = BATCH_USE,
        rescale = config["rescale_for_batch_correction"], # condition
        mnn_fast = config["mnn_use_fast"], # condition
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        hvgs_for_batch_correction = config["hvgs_for_batch_correction"],# condition
        nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], # number
        seeds_umap = config["values"]["batch_correction"]["seeds_umap"]
    script:
        "scripts/01_mnncorrect.R"
        
"""
batch correct using Seurat3

Seurat3 is good at batch correction and among best for multiple batch i
ntegration but not great at recovering DEGs from batch corrected data
unbalanced towards stronger batch effect removal, but successful at 
removing species batch effects. Requires shared cell types between batches 
but no labels, scaling little effect. 
It is easiest to use HVGs and Renormalize using seurat functions.
"""
if BATCH_USE != "Object_ID":
  rule run_seurat3:
      input:
        sce_input = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/fractions/sce_{fraction}-08"
      output:
          sce_output = OUTPUT_DAT + "01_srt3/sce_{fraction}_" + BATCH_USE + "_" + SEURAT_RED + "-01"
      params:
          batch_use = BATCH_USE,
          nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], # number
          nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
          seurat_reduction = SEURAT_RED
      script:
          "scripts/02_seurat.R"
          
#-------------------------------------------------------------------------------

"""
Make batch correction reports 
"""

if BATCH_USE == "Object_ID":
  sce_report_input = rules.run_mnncorrect.output,
if BATCH_USE != "Object_ID": # seurat cannot deal with Object_ID as batch
  sce_report_input = [rules.run_mnncorrect.output, rules.run_seurat3.output]

print(BATCH_USE)

# takes >4h for hscs
rule make_reports:
    input:
        sce_input_raw = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/fractions/sce_{fraction}-08",
        sce_input_corrected = sce_report_input
    output:
        OUTPUT_REP + "batch_correction_report_{fraction}_" + BATCH_USE + ".html"
    params:
        batch_use = BATCH_USE,
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        color_tables = TABLES_PATH
    script:
        "batch_correction_report.Rmd" 
   
print(expand(rules.run_mnncorrect.output, fraction = fractions))
rule make_summary_report:
    input:
        sce_02 = expand(rules.run_mnncorrect.output, fraction = fractions)
    params:
        fractions = fractions
    output:
        OUTPUT_REP + "batch_correction_summary_" + BATCH_USE + ".html"
    script:
        "batch_correction_summary.Rmd" 
