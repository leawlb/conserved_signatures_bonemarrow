#!/bin/python 

"""
Info in batch correction methods from:
Luecken et al. "Benchmarking atlas-level data integration in single-cell genomics", Nat Met 2022
Tran, Ang, Chevrier, Zhang et al. "A benchmark of batch effect correction methods for songle-cell RNA sequencing data", Genome Biologz 2020
"""
#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
METADATA_PATH = config["metadata"]["raw"]
OUTPUT_BASE_PATH = config["paths"]["output_dir"]

# objects from config
IDENTIFIERS = config["metadata"]["identifiers"]
METADATA = pd.read_csv(config["metadata"]["raw"])
VALUES =  config["metadata"]["values"]
BATCH_USE = config["metadata"]["batch_use"] # which SCE col to use as batch

#-------------------------------------------------------------------------------

targets = METADATA["Species_ID"]
targets = targets.drop_duplicates()
targets = targets.squeeze()
targets = targets.tolist()
targets = targets + "all" # for the SCE object containing all samples
print(targets)

targets_06 = [OUTPUT_BASE_PATH + "/06_mrge/" + "sce_" + x + "-06" for x in targets]
targets_07 = [OUTPUT_BASE_PATH + "/07_rnrm/" + "sce_" + x + "-07" for x in targets]
# construct specific targets for each used batch and each used method only if method is to be used
targets_08a = [OUTPUT_BASE_PATH + "/08_mnncorrect/" + "sce_" + x + "_correctedby_" + BATCH_USE + "-08" for x in targets] if run_mnncorrect else []
targets_08b = [OUTPUT_BASE_PATH + "/08_seurat3/" + "sce_" + x + "_correctedby_" + BATCH_USE + "-08" for x in targets] if run_seurat3 else []
targets_08c = [OUTPUT_BASE_PATH + "/08_scmerge/" + "sce_" + x + "_correctedby_" + BATCH_USE + "-08" for x in targets] if run_scmerge else []

targets_08 = targets_08a + targets_08b + targets_08c

#-------------------------------------------------------------------------------

localrules: all  

rule_all:
    input:
        targets_06,
        targets_07,
        targets_08

"""
Merge all datasets into one big dataset, and four species-specific SCE objects
Input data is located in 02_preprocessing/04_norm/
"""
rule merge_datasets:
    input:
        sce_04_path = config["paths"]["preprocessing"] + "/norm/", 
    output:
        sce_06_path = OUTPUT_BASE_PATH + "/06_mrge/" # rest is specified in script
    params:
        targets = targets
    script:
        "scripts/06_merge_datasets.R" 

"""
# quick and basic renormalization, scaling, HVG calculation

Scaling improves batch effect removal but worsens bioconservation
HVG selection improves performance but restricts analysis
"""
# quick and basic renormalization, scaling, HVG calculation
rule renormalize:
    input:
        sce_06 = rule.merge_data.output + "sce_{species}-06"
    output:
        sce_07 = OUTPUT_BASE_PATH + "/rnrm/sce_{species}-07"
        hvgs = OUTPUT_BASE_PATH + "/rnrm/hvg_{species}-07"
    params:
        batch_use = BATCH_USE
    script:
        "scripts/07_renormalize.R"

"""
batch correct using MNNcorrect

MNNcorrect is good at recovering DEGs from batch corrected data and bioconservation,
but slow and does not perform well on batch correction 
Requires shared cell types between batches but no labels
(fastMNN) seems to balance vatch effect removal and bioconservation
"""
if config["run_mnncorrect"]:
    rule renormalize:
        input:
            sce_07 = rule.renormalize.output 
        output:
            sce_08 = OUTPUT_BASE_PATH + "/08_mnncorrect/sce_{species}_correctedby_" + BATCH_USE + "-08"
        params:
            hvgs_for_batch_correction = config["hvgs_for_batch_correction"]
        script:
            "scripts/08_mnncorrect.R"

"""
batch correct using Seurat3

Seurat3 is good at batch corection and among best for multiple batch integration
but not great at recovering DEGs from batch corrected data
unbalanced towards stronger batch effect removal, but successful at removing species batch effects
Requires shared cell types between batches but no labels, scaling little effect
"""
if config["run_seurat3"]:
    rule renormalize:
        input:
            sce_07 = rule.renormalize.output 
        output:
            sce_08 = OUTPUT_BASE_PATH + "/08_seurat3/sce_{species}_correctedby_" + BATCH_USE + "-08"
        params:
            hvgs_for_batch_correction = config["hvgs_for_batch_correction"]
         script:
            "scripts/08_seurat3.R"

"""
batch correct using scMerge

scMerge and among best for multiple batch integration,
is ok but not great at batch corection and recovering DEGs
seems balanced but is also slow
"""
if config["run_scmerge"]:
    rule renormalize:
        input:
            sce_07 = rule.renormalize.output 
        output:
            sce_08 = OUTPUT_BASE_PATH + "/08_scmerge/sce_{species}_correctedby_" + BATCH_USE + "-08"
        params:
            hvgs_for_batch_correction = config["hvgs_for_batch_correction"]
         script:
            "scripts/08_scmerge.R"
