#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
OUTPUT_BASE = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]
VALUES = config["values"]["clustering"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno/"
OUTPUT_REP = OUTPUT_BASE + "/reports/02_sce_anno/02_clustering/"

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

# construct paths for all possible outputs/targets, required for rule all
targets = []
for f in fractions:
  #targets = targets + [OUTPUT_DAT + "02_clst/hrchl_clust/sce_" + f + "-02"]
  #targets = targets + [OUTPUT_DAT + "02_clst/seurt_clust/sce_" + f + "-02"]
  targets = targets + [OUTPUT_DAT + "02_clst/louvn_clust/sce_" + f + "-02"]
  targets = targets + [OUTPUT_REP + "clustering_report_full_" + f +".html"]
  if config["calculate_k"]:
    targets = targets + [OUTPUT_REP + "clustering_report_k_" + f +".html"]
    
#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------

rule hierarchical_clustering:
    input: 
        sce_input = OUTPUT_DAT + "01_mnnc/sce_{fraction}_Object_ID-01"
    output:
        sce_output = OUTPUT_DAT + "02_clst/hrchl_clust/sce_{fraction}-02"
    params:
        number_k = VALUES["number_k"],
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        size_subs_hscs = VALUES["size_subs_hscs"]
    script:
        "scripts/02_hierarchical_clustering.R"
        
rule seurat_clustering:
    input: 
        sce_input = OUTPUT_DAT + "01_mnnc/sce_{fraction}_Object_ID-01"
    output:
        sce_output = OUTPUT_DAT + "02_clst/seurt_clust/sce_{fraction}-02"
    params:
        resolution = VALUES["resolution"]
    script:
        "scripts/02_seurat_clustering.R"
        
rule louvain_clustering:
    input: 
        sce_input = OUTPUT_DAT + "01_mnnc/sce_{fraction}_Object_ID-01"
    output:
        sce_output = OUTPUT_DAT + "02_clst/louvn_clust/sce_{fraction}-02"
    params:
        k_louvain = VALUES["k_louvain"]
    script:
        "scripts/02_louvain_clustering.R"

#-------------------------------------------------------------------------------   

if config["calculate_k"]:     
  rule make_report_k:
      input:
          sce_input = OUTPUT_DAT + "01_mnnc/sce_{fraction}_Object_ID-01"
      output:
          OUTPUT_REP + "clustering_report_k_{fraction}.html"
      params: 
          size_subs_hscs = VALUES["size_subs_hscs"],
      threads:
          20
      script:
          "clustering_report_k.Rmd"

rule make_report_full:
    input:
        sce_input = OUTPUT_DAT + "01_mnnc/sce_{fraction}_Object_ID-01",
        #sce_h = rules.hierarchical_clustering.output,
        #sce_s = rules.seurat_clustering.output,
        sce_l = rules.louvain_clustering.output
    params:
        color_tables = TABLES_PATH,
        size_subs_hscs = VALUES["size_subs_hscs"],
        number_k = VALUES["number_k"]
    output:
        OUTPUT_REP + "clustering_report_full_{fraction}.html"
    script:
        "clustering_report_full.Rmd"
