#!/bin/python 

# Prepping SCEs for CCI calculation

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/02_cci_calc/"
OUTPUT_REP = OUTPUT_BASE + "/reports/04_cci_calc/"

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
    targets = targets + [OUTPUT_DAT + "01_perc/perc_df_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "02_intm/interaction_matrix_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "02_intm/datasheet_" + s + "_" + a]


#targets = targets + [OUTPUT_REP + "sce_downsampling_report.html"]
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
 
#-------------------------------------------------------------------------------


""" 
At this point, certain conditions in the SCE objects must be met.
These conditions are fulfilled in 01_cci_prep.

- cluster or cell type labels:
  - should be in colData slot called "Identity"
  - should not be able to grep each other (e.g. "B cell" and "pre-B cell")
  - should not contain ungreppable characters (+, (, ), or &)
  - should be factored in a preferable sequence

- rownames should contain Gene Symbols
- a rowData slot should countain Gene IDs in slot "ENSMUS_ID"

- Emitter/Receiver Assignment:
  - should be in colData slot called "Assignment"
  - should contain "emitter" for designated emitter cells
  - should contain "receiver" for designated receiver cells
  
- Quality
  - all SCE objects to be compared should have comparable counts matrices
  - downsampled to lowest quality samples

Please refer to 01_cci_prep for detailed steps to ensure these conditions.

"""

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
Extract percentage of cells per identity in which each gene is expressed

This is important because genes expressed in < X percent increase noise.

Output = dataframe containing % of expressing cells for each gene per identity

"""

rule get_perc_df:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_prep/05_down/sce_{species}_{age}-05"
    output:
        perc_df = OUTPUT_DAT + "01_perc/perc_df_{species}_{age}"
    script:
        "scripts/01_get_perc_df.R" 

#-------------------------------------------------------------------------------

"""
Extract interaction matrix

This is an important step in the construction of CCI objects.
Only genes that are expressed in at least X% of cells per identity are
considered for CCI calculation.

From the counts matrix of SCE objects, only ligand gene counts are extracted for 
emitters, and only receptor genes for receivers.

Output = An interaction matrix with cols = cells and rows = interactions

"""

rule get_interaction_mat:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_prep/05_down/sce_{species}_{age}-05",
        perc_df = rules.get_perc_df.output,
        lrdb = config["metadata"]["path"] + "/CCI/lrdb"
    output:
        interaction_mat = OUTPUT_DAT + "02_intm/interaction_matrix_{species}_{age}"
        datasheet = OUTPUT_DAT + "02_intm/datasheet_{species}_{age}"
    params:
        min_perc = config["values"]["cci_calc"]["min_perc"]
    script:
        "scripts/02_get_interaction_mat.R" 
        
#-------------------------------------------------------------------------------

"""
Interaction Ranking
"""
