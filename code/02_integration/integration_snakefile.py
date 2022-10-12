#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
METADATA_PATH = config["metadata"]["raw"]
OUTPUT_BASE_PATH = config["paths"]["integration"]

# objects from config
IDENTIFIERS = config["metadata"]["identifiers"]
METADATA = pd.read_csv(config["metadata"]["raw"])
VALUES =  config["metadata"]["values"]

#-------------------------------------------------------------------------------

rule_all

"""
Merge all datasets into one big, and four species-specific SCE objects
Input data is located in 01_preprocessing/norm
"""
rule merge_datasets:
    input:
        sce_04_path = config["paths"]["preprocessing"] + "/norm/", 
    output:
        sce_06_path = OUTPUT_BASE_PATH + "/mrge/" 
        # rest of path is specified in 01_integration.R
    params:
        targets = targets
        run_sce_all = config["run_sce_all"] # specify if big SCE should be made
    script:
        "scripts/01_integration.R" 

rule renormalize:
    input:
        sce_06 = rule.merge_data.output + "sce_{species}-06"
    output:
        sce_07 = OUTPUT_BASE_PATH + "/rnrm/sce_{species}-07"
    params:
        renorm_method = config["methods"]["renorm_method"]
    script:
        "scripts/02_renormalize.R"

