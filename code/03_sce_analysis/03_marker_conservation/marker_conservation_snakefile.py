#!/bin/python 

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
print(OUTPUT_BASE)

RUN_AGE_COMP = config["run_age_comparison"]

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_str.RData"]
targets = targets + [OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_hsc.RData"]

if RUN_AGE_COMP:
  targets = targets + [OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_str.RData"]
  targets = targets + [OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_hsc.RData"]
  targets = targets + [OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/03_markers/plot_marker_age_comparison.html"]

#-------------------------------------------------------------------------------

localrules: all  

rule all: 
  input:
      targets

#-------------------------------------------------------------------------------

rule marker_per_sp_per_celltype:
    input:
        data_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10",
        data_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10"
    resources:
        mem_mb=50000,
          queue = "long-debian"
    threads: 6
    output:
        output_str = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_str.RData",
        output_hsc = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_hsc.RData"
    script:
        "01_marker_per_sp_per_celltype.R"

# signatures are only calculated a step later
# so run this rule only when 04 was already done

if RUN_AGE_COMP:

  rule marker_per_age_comparison:
      input:
          data_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
          data_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10",
          markers_cons_hsc = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own/01_gens/geneset_list_hsc",
          markers_cons_str = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own/01_gens/geneset_list_str"
      resources:
          mem_mb=50000,
          queue = "long-debian"
      threads: 30
      output:
          output_str = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_str.RData",
          output_hsc = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_hsc.RData"
      script:
          "02_marker_age_comparison.R"
          
          
  rule report:
      input:
          markers_cons_hsc =  OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_hsc.RData",
          markers_cons_str =  OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_str.RData",
          markers_conservation_hsc =  OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_hsc.RData",
          markers_conservation_str =  OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_str.RData"
      output: 
          OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/03_markers/plot_marker_age_comparison.html"
      resources:
          mem_mb=30000,
          queue = "short-debian"
      threads: 2
      script:
          "plot_marker_age_comparison.Rmd"
