#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
OUTPUT_BASE = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno/"
OUTPUT_REP = OUTPUT_BASE + "/reports/02_sce_anno/03_anno_clusters/"


#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")

# cannot take from metadata since cluster number is dependent on fraction
clusters_hsc = list(range(1,14))
clusters_hsc = list(map(str,clusters_hsc))

clusters_str = list(range(1,18))
clusters_str = list(map(str,clusters_str))

# construct paths for all possible outputs/targets, required for rule all
targets = []
for c in clusters_hsc:
  targets = targets + [OUTPUT_DAT + "04_sepd/sce_hsc_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_REP + "annotation_hsc_cluster_" + c + ".html" ] 

for c in clusters_str:
  targets = targets + [OUTPUT_DAT + "04_sepd/sce_str_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_REP + "annotation_str_cluster_" + c + ".html"] 

for f in fractions:
  targets = targets + [OUTPUT_DAT + "05_annc/01_markers/markers_" + f]
  targets = targets + [OUTPUT_DAT + "05_annc/02_goan/go_" + f]
  #targets = targets + [OUTPUT_DAT + "05_annc/03_sce/sce_" + f + "-05"]

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
    cluster="[0-9]+",

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets

#-------------------------------------------------------------------------------
# dummy rule to allow using clusters as wildcard
output = expand(OUTPUT_DAT + "04_sepd/sce_hsc_cluster_{cluster}-sep", cluster = clusters_hsc)
output = output + expand(OUTPUT_DAT + "04_sepd/sce_str_cluster_{cluster}-sep", cluster = clusters_str)
rule separate_sce:
    input: 
        sce_input = expand(OUTPUT_DAT + "03_clst/louvn_clust/sce_{fraction}-03", fraction = fractions)
    output:
        sce_output = output
    script:
        "scripts/04_separate_dummy.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
annotation of clusters
"""

# get marker genes for each cluster and do GO analysis on the marker genes.
rule find_markers:
    input: 
        sce_input = OUTPUT_DAT + "03_clst/louvn_clust/sce_{fraction}-03"
    output:
        markers = OUTPUT_DAT + "05_annc/01_markers/markers_{fraction}"
    params:
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
    script:
        "scripts/05_markers_clusters.R"

# perform preliminary GO analysis for overview (Very basic)
rule go_analysis:
    input: 
        markers = rules.find_markers.output
    output:
        go = OUTPUT_DAT + "05_annc/02_goan/go_{fraction}"
    script:
        "scripts/05_go_clusters.R"

# report on marker gene expression and GO 
rule make_report:
    input: 
        markers = rules.find_markers.output,
        go = rules.go_analysis.output,
        sce_sep = OUTPUT_DAT + "04_sepd/sce_{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "03_clst/louvn_clust/sce_{fraction}-03",
        gene_list_all = config["metadata"]["path"] + "/gene_list_{fraction}_all.csv"
    output:
        OUTPUT_REP + "annotation_{fraction}_cluster_{cluster}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "anno_clusters_report.Rmd"
        
# assign annotation labels informed from report
rule assign_annotation:
    input:
        sce_input = OUTPUT_DAT + "03_clst/louvn_clust/sce_{fraction}-03",
        #anno_list =  config["metadata"]["path"] + "/anno_list_{fraction}.csv"
    output:      
        sce_output = OUTPUT_DAT + "05_annc/03_sce/sce_{fraction}-05"
    script:
        "scripts/05_anno_clusters.R"
      
