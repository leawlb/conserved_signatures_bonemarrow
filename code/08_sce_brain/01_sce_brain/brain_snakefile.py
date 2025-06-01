#!/bin/python 

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + "/data/scRNAseq"
print(OUTPUT_BASE)

INPUT_DATASET = config["base"] + config["metadata_paths"]["dataset_brain_path"]
print(INPUT_DATASET)

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/01_list_nDEGs_all.rds"]
targets = targets + [OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/02_marker_expression_primates.rds"]
targets = targets + [OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/02_marker_conserved_primates.rds"]
targets = targets + [OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/03_resolution_scores_markers.rds"]
targets = targets + [OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/03_resolution_scores_signature.rds"]
targets = targets + [OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/04_recluster/mouse_reclust_hs.rds"]

targets = targets + [OUTPUT_BASE + "/main_analysis/sce_objects/reports/08_sce_brain/report.html"]

print(targets)

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
  input:
      targets

#-------------------------------------------------------------------------------

rule nDEGs:
    input:
        data_input = INPUT_DATASET
    resources:
        mem_mb = 20000,
        queue = "medium-debian"
    threads: 15
    output:
        data_output = OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/01_list_nDEGs_all.rds"
    script:
        "01_nDEGs.R"

rule markers:
    input:
        data_input = INPUT_DATASET
    resources:
        mem_mb = 20000,
        queue = "medium-debian"
    threads: 15
    output:
        markers_conservation = OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/02_marker_expression_primates.rds",
        cons_markers = OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/02_marker_conserved_primates.rds"
    script:
        "02_markers.R"

rule resolution:
    input:
        data_input = INPUT_DATASET,
        cons_markers = rules.markers.output.cons_markers,
        nDEGS = rules.nDEGs.output.data_output
    resources:
        mem_mb = 90000,
        queue = "medium-debian"
    threads: 15
    output:
        all_scores_markers = OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/03_resolution_scores_markers.rds",
        all_scores_signature = OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/03_resolution_scores_signature.rds",
    script:
        "03_resolution.R"
 
rule recluster:
    input:
        data_input = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS",
        cons_markers = rules.markers.output.cons_markers,
        nDEGS = rules.nDEGs.output.data_output,
        markers_conservation = rules.markers.output.markers_conservation,
    params:
        base_path = OUTPUT_BASE,
        brain_path = "/main_analysis/sce_objects/08_sce_brain/"
    resources:
        mem_mb = 40000,
        queue = "medium-debian"
    threads: 15
    output:
        species_human_only = OUTPUT_BASE + "/main_analysis/sce_objects/08_sce_brain/04_recluster/mouse_reclust_hs.rds"
    script:
        "04_recluster.R"


"""
rule quantify_reclust

# 05_quantify_reclust is mostly loading and plotting, so no output 
aside from a pdf.

Copied in its entirety into the report.

"""

rule report:
    input:
        all_scores_markers = rules.resolution.output.all_scores_markers,
        all_scores_signature = rules.resolution.output.all_scores_signature,
        species = rules.recluster.output.species_human_only,
        cons_markers = rules.markers.output.cons_markers,
        nDEGS = rules.nDEGs.output.data_output
    params:
        base_path = OUTPUT_BASE,
        brain_path = "/main_analysis/sce_objects/08_sce_brain/"
    resources:
        mem_mb = 20000,
        queue = "medium-debian"
    threads: 2
    output:
        OUTPUT_BASE + "/main_analysis/sce_objects/reports/08_sce_brain/report.html"
    script:
        "report.Rmd"
        
