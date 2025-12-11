#!/bin/python 

"""
This script is to perform signature-based re-clustering followed by blind
annotation, but this time with batch correction.

The chosen datasets are ts_bone_marrow and ts_hscs_progenitors, because
they are both shown in figure 4 (as of now) and are good examples of data with 
strong technical or batch effects due to use of different sex and assays.

Because it is much simpler and each dataset needs to be treated 
individually, I will perform the analysis in .Rmd scripts in contrast to
the usual .R scripts, one for each dataset.
Then, I will save the object from within the .Rmd script.
I am also loading the required datasets within each separate .Rmd script.

See also the analysis parts 03_custom and 04_blind_anno where I use the
same concept.

"""

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/06_batchcorrection"
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures"

COLORS = config["base"] + config["metadata_paths"]["colors"]

METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_REP + "/bc_ts_bonemarrow.html"]


#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# tabula sapiens adult bone marrow
# using the same environment as for the original signature-based
# clustering, which is the snkmk_isbm env

rule batchcorrection_ts_bonemarrow:
    input: 
        seu_input = config["base"] + config["metadata_paths"]["datasets_other_path"] + "/ts_bone_marrow",
        ensembl_sign = OUTPUT_DAT + "/01_reclustering_own/02_endf/ensembl_sign_hsc",
        resolution_df_path = config["base"] + config["metadata_paths"]["resolution_other"]
    output:
        OUTPUT_REP + "/bc_ts_bonemarrow.html"
    resources:
        mem_mb = 50000,
        queue = "medium-debian"
    threads: 10
    params:
        plotting = "../../source/plotting.R",
        OUTPUT_DAT = OUTPUT_DAT
    script:
        "06_include_batchcorrection/bc_ts_bonemarrow.Rmd"

