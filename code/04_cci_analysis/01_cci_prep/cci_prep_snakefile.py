#!/bin/python 

# Prepping SCEs for CCI calculation

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_prep/"
OUTPUT_REP = OUTPUT_BASE + "/reports/04_cci_prep/"

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")

targets = []

for s in species:
  for a in age:
    targets = targets + [OUTPUT_DAT + "01_sprp/sce_" + s + "_" + a + "-01"]
    targets = targets + [OUTPUT_DAT + "02_asnm/sce_" + s + "_" + a + "-02"]
    targets = targets + [OUTPUT_DAT + "03_down/sce_" + s + "_" + a + "-03"]

targets = targets + [OUTPUT_REP + "sce_downsampling_report.html"]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
 
#-------------------------------------------------------------------------------


""" 
PREP SCE

The SCE objects need to be prepared before CCI calculation
"""

# add cluster labels from annotated fraction merges to age merges
rule add_clusterlabels:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/ages/sce_{species}_{age}-08",
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/02_clst/louvn_clust/sce_hsc-test",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/02_clst/louvn_clust/sce_str-test"
    output:
        sce_output = OUTPUT_DAT + "01_sprp/sce_{species}_{age}-01"
    script:
        "scripts/01_add_clusterlabels.R" 

#-------------------------------------------------------------------------------

# add assignments for each cell type (emitter, receiver or remove)
# remove cells or cell types that are not wanted or required for CCI calculation
# or that are not abundant enough
rule celltype_assignment:
    input:
        sce_input = rules.add_clusterlabels.output,
        assignment = "assignment.txt"
    params:
        min_cells = config["values"]["cci_prep"]["min_cells"]
    output:
        sce_output = OUTPUT_DAT + "02_asnm/sce_{species}_{age}-02"
    script:
        "scripts/02_add_assignment.R" 
        
# downsample to the lowest sample for comparability between samples
rule downsample:
    input:
        sce_input_path = expand(rules.celltype_assignment.output, species = species, age = age)
    output:
        sce_output_path = expand(OUTPUT_DAT + "03_down/sce_{species}_{age}-03", species = species, age = age)
    script:
        "scripts/03_downsample.R" 

# sce prep report
rule make_downsampling_reports:
    input:
        sce_input_path =  expand(rules.downsample.output, species = species, age = age)
    output:
        OUTPUT_REP + "sce_downsampling_report.html"
    script:
        "sce_downsampling_report.Rmd"
        
#-------------------------------------------------------------------------------

# prepare LRDB
# CellTalkDB mouse_lr_pair.rds downloaded from https://github.com/ZJUFanLab/CellTalkDB/blob/master/database/mouse_lr_pair.rds
# on Jul 6, 2023 at 2:50PM






## PREP LRDB


