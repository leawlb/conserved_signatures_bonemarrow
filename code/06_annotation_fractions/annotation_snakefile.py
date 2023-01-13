#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import sys 

# paths from config
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["raw"])
print(OUTPUT_BASE_PATH)

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
print(species) 

# construct paths for all possible outputs/targets, required for rule all
targets = [OUTPUT_BASE_PATH + "/sce_objects/14_scmp_cluster/sce_" + s + "-14" for s in species]

for s in species:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/14_scmp_cluster/sce_" + s + "_nbc-14"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/15_scmp_reference/sce_" + s + "-15"]
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/06_annotation/" + s + "/annotation_species_report_" + s + "_refanno.html"]

if config["run_annotation_summary_scmap"]:
  targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/reports/06_annotation/annotation_summary_scmap.html"]

print(targets)

#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------


rule cluster_scmap_annotation:
    input: 
        sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge/sce_{species}-07",
        sce_12 = OUTPUT_BASE_PATH + "/sce_objects/12_lcls/sce_{species}-12",
        mmus = OUTPUT_BASE_PATH + "/sce_objects/12_lcls/sce_mmus-12"
    output:
        sce_14 = OUTPUT_BASE_PATH + "/sce_objects/14_scmp_cluster/sce_{species}-14",
        sce_14_nbc = OUTPUT_BASE_PATH + "/sce_objects/14_scmp_cluster/sce_{species}_nbc-14"
    params:
        species = species
    script:
        "scripts/14_scmap_mmus_cluster_annotation.R"

rule reference_scmap_annotation:
    input: 
        sce_14 = rules.cluster_scmap_annotation.output.sce_14
    output:
        sce_15 = OUTPUT_BASE_PATH + "/sce_objects/15_scmp_reference/sce_{species}-15"
    params:
        ref_baccin_sce = config["metadata"]["ref_baccin_sce"],
        ref_dahlin_sce = config["metadata"]["ref_dahlin_sce"],
        ref_dolgalev_sce = config["metadata"]["ref_dolgalev_sce"]
    script:
        "scripts/15_scmap_reference_celltype_annotation.R"
        
rule make_species_report:
    input: 
        sce_15 = rules.reference_scmap_annotation.output
    output:
        OUTPUT_BASE_PATH + "/sce_objects/reports/06_annotation/{species}/annotation_species_report_{species}_refanno.html"
    params:
         color_tables = TABLES_PATH
    script:
        "annotation_species_report_refanno.Rmd"
 
if config["run_annotation_summary_scmap"]:
  rule make_summary_report_scmap:
      input:
          sce_14_path = OUTPUT_BASE_PATH + "/sce_objects/14_scmp_cluster/"
      params:
         color_tables = TABLES_PATH,
         species = species
      output:
          OUTPUT_BASE_PATH + "/sce_objects/reports/06_annotation/annotation_summary_scmap.html"
      script:
          "annotation_summary_scmap.Rmd" 
