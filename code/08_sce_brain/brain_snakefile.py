#!/bin/python 

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + "/data"
print(OUTPUT_BASE)

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/01_list_nDEGs_all.rds"]
targets = targets + [OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/02_marker_expression_primates.rds"]
targets = targets + [OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/02_marker_conserved_primates.rds"]
targets = targets + [OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/03_resolution_scores.rds"]

print(targets)

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
  input:
      targets

#-------------------------------------------------------------------------------

rule nDEGs:
    input:
        data_input = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS"
    resources:
        mem_mb = 20000,
        queue = "medium"
    threads: 15
    output:
        data_output = OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/01_list_nDEGs_all.rds"
    script:
        "01_nDEGs.R"

rule markers:
    input:
        data_input = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS"
    resources:
        mem_mb = 20000,
        queue = "medium"
    threads: 15
    output:
        markers_conservation = OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/02_marker_expression_primates.rds",
        core_markers = OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/02_marker_conserved_primates.rds"
    script:
        "02_markers.R"

rule resolution:
    input:
        data_input = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS",
        core_markers = rules.markers.output.core_markers
    resources:
        mem_mb = 20000,
        queue = "medium"
    threads: 15
    output:
        all_scores = OUTPUT_BASE + "/scRNAseq/main_analysis/sce_objects/08_sce_brain/03_resolution_scores.rds",
    script:
        "03_resolution.R"
