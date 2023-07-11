#!/bin/python 

# Prepping SCEs for CCI calculation

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/01_cci_preparation/"

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
    targets = targets + [OUTPUT_DAT + "06_perc/perc_df_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "07_intm/interaction_matrix_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "07_intm/datasheet_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "08_rank/interaction_ranking_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "09_intl/interaction_list_" + s + "_" + a]

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
  - should not be able to grep each other (e.g. "BC" is containted in "pre-BC")
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

This is important because genes expressed in <X percent of cells increase noise.

Output = dataframe containing % of expressing cells for each gene per identity

"""

rule get_perc_df:
    input:
        sce_input = OUTPUT_DAT + "05_down/sce_{species}_{age}-05"
    output:
        perc_df = OUTPUT_DAT + "06_perc/perc_df_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/06_get_perc_df.R" 

#-------------------------------------------------------------------------------

"""
Extract interaction matrix

This is an important step in the construction of CCI objects.
Only genes that are expressed in at least X% of cells per identity are
considered for CCI calculation.

From the counts matrix of SCE objects, only ligand gene counts are extracted for 
emitters, and only receptor genes for receivers.

Main output = An interaction matrix with cols = cells and rows = interactions

"""

rule get_interaction_mat:
    input:
        sce_input = OUTPUT_DAT + "05_down/sce_{species}_{age}-05",
        perc_df = rules.get_perc_df.output,
        lrdb = config["metadata"]["path"] + "/CCI/lrdb"
    output:
        interaction_mat = OUTPUT_DAT + "07_intm/interaction_matrix_{species}_{age}",
        datasheet = OUTPUT_DAT + "07_intm/datasheet_{species}_{age}"
    params:
        min_perc = config["values"]["cci_calc"]["min_perc"],
        main_functions = "../../source/cci_functions_prep_main.R",
        help_functions = "../../source/cci_functions_prep_help.R"
    script:
        "scripts/07_get_interaction_mat.R" 
        
#-------------------------------------------------------------------------------

"""
Interaction Ranking

The heart of this method.

- For each emitter, mean ligand gene expression across cells is ranked 
- For each receiver, mean receptor gene expression across cells is ranked. 

- For each emitter-receiver cell type pair, ligand-receptor interactions are
  then scored. Ligand and receptor ranks are summed, then normalized.

- To reduce the influence of the number of detected ligands and receptors per
  cell type on the ranking, ranks are levelled to a "top_level". 
  As such, the top level in each cell type will always be X. 

Main output = List of matrices, including scores, ligand rank, and receptor rank 
information.

"""

rule interaction_ranking:
    input:
        interaction_mat = rules.get_interaction_mat.output.interaction_mat,
        datasheet = rules.get_interaction_mat.output.datasheet
    output:
        interaction_ranking = OUTPUT_DAT + "08_rank/interaction_ranking_{species}_{age}"
    params:
        top_level = config["values"]["cci_calc"]["top_level"],
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/08_interaction_ranking.R" 
        
"""

Making the interaction list.

The ranking and scoring are transformed into a more convenient format.

Output = list of annotated objects with all relevant information and all
data required for CCI analysis.
"""

rule interaction_list:
    input:
        interaction_ranking = rules.interaction_ranking.output,
        lrdb = config["metadata"]["path"] + "/CCI/lrdb"
    output:
        interaction_list = OUTPUT_DAT + "09_intl/interaction_list_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_prep_main.R"
    script:
        "scripts/09_get_interaction_list.R" 
        

## QUALITY CONTROL will be performed in 03_metrics
