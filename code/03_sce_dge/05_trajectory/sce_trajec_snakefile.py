#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/05_trajectory"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/05_trajectory"

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
COLORS = config["base"] + config["metadata_paths"]["colors"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
print(METADATA)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")
fractions = ["hsc"]

#-------------------------------------------------------------------------------

# construct paths for all possible outputs/targets, required for rule all
targets = []

for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_mprp/cds_" + f]
  #targets = targets + [OUTPUT_DAT + "/11_mprp/pdf1_" + f + ".pdf"]
  #targets = targets + [OUTPUT_DAT + "/11_mprp/pdf2_" + f + ".pdf"]
  targets = targets + [OUTPUT_DAT + "/02_mgrp/cds_" + f]
  #targets = targets + [OUTPUT_DAT + "/12_mgrp/pdf_" + f + ".pdf"]
  targets = targets + [OUTPUT_DAT + "/03_pstm/cds_" + f]
#  for s in species:
#    targets = targets + [OUTPUT_DAT + "/10_ssce/sce_" + f + "_" + s]
#    targets = targets + [OUTPUT_DAT + "/11_mprp/cds_" + f + "_" + s]
#    targets = targets + [OUTPUT_DAT + "/12_mgrp/cds_" + f + "_" + s]
#    targets = targets + [OUTPUT_DAT + "/12_mgrp/pdf_" + f + "_" + s + ".pdf"]
    #targets = targets + [OUTPUT_DAT + "/13_pstm/cds_" + f + "_" + s]

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
    fraction="[a-z]+",

localrules: all  

rule all: 
    input:
        targets
        
        
#-------------------------------------------------------------------------------
"""
# Convert SCE to CDS (monocle3) objects and use default monocle3 pipeline
# to calculate trajectories for both fractions.

Pipeline from: https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

On all HSPCs, removing cell cycle genes
"""

# all steps before finding the graph, which is required for finding the nodes
rule monocle3_prep:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        cell_cycle_genes = config["base"] + config["metadata_paths"]["cell_cycle_genes"]
    output:
        cds_output = OUTPUT_DAT + "/01_mprp/cds_{fraction}",
        #pdf_output1 = OUTPUT_DAT + "/11_mprp/pdf1_{fraction}.pdf",
        #pdf_output2 = OUTPUT_DAT + "/11_mprp/pdf2_{fraction}.pdf"
    params:
        k_louvain_monoc = config["values"]["02_sce_anno"]["k_louvain"]
    conda: 
        "monocle3.yaml" 
    script:
        "scripts/01_monocle3.R"
        
# find the graph and export nodes to be chosen      
rule find_graph:
    input:
        cds_input = rules.monocle3_prep.output.cds_output
    output:
        cds_output = OUTPUT_DAT + "/02_mgrp/cds_{fraction}",
        #pdf_output = OUTPUT_DAT + "/12_mgrp/pdf_{fraction}.pdf"
    conda: 
        "monocle3.yaml" 
    script:
        "scripts/02_find_graph.R"

# order cell by pseudotime  
rule pseudotime:
    input:
        cds_input = rules.find_graph.output.cds_output
    output:
        cds_output = OUTPUT_DAT + "/03_pstm/cds_{fraction}"
    conda: 
        "monocle3.yaml" 
    script:
        "scripts/03_pseudotime.R"
        
