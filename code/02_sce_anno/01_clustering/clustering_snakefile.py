#!/bin/python 

"""
Batch correction, clustering and cluster annotation

This is the version for tested pipelines. 
For testing pipelines with other methods see 02_sce_anno_TEST

"""
#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
OUTPUT_BASE = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
METADATA = pd.read_csv(config["metadata"]["table"])
BATCH_USE = config["values"]["batch_correction"]["batch_use_fractions"] # which colData to use as batch
# seurat kept for QC sake but not used for downstream analysis
SEURAT_RED = config["values"]["batch_correction"]["seurat_reduction"] # which seurat method to use

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno/"
OUTPUT_REP = OUTPUT_BASE + "/reports/02_sce_anno/"

# cannot take from metadata since cluster number is dependent on fraction
clusters_hsc = list(range(1,13))
clusters_hsc = list(map(str,clusters_hsc))

clusters_str = list(range(1,18))
clusters_str = list(map(str,clusters_str))

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
fractions = get_list(metadata = METADATA, column = "Fraction_ID")

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "01_mnnc/sce_" + f + "_" + BATCH_USE + "-01"]
  targets = targets + [OUTPUT_DAT + "02_clst/louvn_clust/sce_" + f + "-02"]
  targets = targets + [OUTPUT_REP + "01_batch_correction/batch_correction_report_" + f + "_" + BATCH_USE + ".html"]
  targets = targets + [OUTPUT_REP + "02_clustering/clustering_report_full_" + f +".html"]
  targets = targets + [OUTPUT_DAT + "04_annc/01_markers/markers_" + f]
  targets = targets + [OUTPUT_DAT + "04_annc/02_goan/go_" + f]
  targets = targets + [OUTPUT_DAT + "04_annc/03_sce/sce_" + f + "-04"]

for c in clusters_hsc:
  targets = targets + [OUTPUT_DAT + "03_sepd/sce_hsc_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_REP + "03_anno_clusters/annotation_hsc_cluster_" + c + ".html" ] 

for c in clusters_str:
  targets = targets + [OUTPUT_DAT + "03_sepd/sce_str_cluster_" + c + "-sep"]
  targets = targets + [OUTPUT_REP + "03_anno_clusters/annotation_str_cluster_" + c + ".html"] 

#-------------------------------------------------------------------------------

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
    input:
        targets
 
#-------------------------------------------------------------------------------
"""
Info on batch correction methods from:
Luecken et al. "Benchmarking atlas-level data integration in single-cell 
genomics", Nat Met 2022
Tran, Ang, Chevrier, Zhang et al. "A benchmark of batch effect correction 
methods for single-cell RNA sequencing data", Genome Biology 2020

Batch correction using MNNcorrect:

"""
rule run_mnncorrect:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/fractions/sce_{fraction}-08"
    output:
        sce_output = OUTPUT_DAT + "01_mnnc/sce_{fraction}_" + BATCH_USE + "-01"
    params:
        batch_use = BATCH_USE,
        nr_hvgs_batch_correction = config["values"]["batch_correction"]["nr_hvgs_batch_correction"], 
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
        seeds_umap = config["values"]["batch_correction"]["seeds_umap"]
    script:
        "scripts/01_mnncorrect.R"
        
          
"""
Make batch correction report
"""

rule make_batchcorrection_reports:
    input:
        sce_input_raw = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/fractions/sce_{fraction}-08",
        sce_input_corrected = rules.run_mnncorrect.output
    output:
        OUTPUT_REP + "01_batch_correction/batch_correction_report_{fraction}_" + BATCH_USE + ".html"
    params:
        batch_use = BATCH_USE,
        color_tables = TABLES_PATH
    script:
        "batch_correction_report.Rmd" 
  
        
#-------------------------------------------------------------------------------
"""
Louvain clustering
"""

rule louvain_clustering:
    input: 
        sce_input = rules.run_mnncorrect.output
    output:
        sce_output = OUTPUT_DAT + "02_clst/louvn_clust/sce_{fraction}-02"
    params:
        k_louvain = config["values"]["clustering"]["k_louvain"]
    script:
        "scripts/02_louvain_clustering.R"
  
"""
Make clustering reports 
"""      
        
rule make_clustering_report:
    input:
        sce_input = rules.run_mnncorrect.output,
        sce_l = rules.louvain_clustering.output
    params:
        color_tables = TABLES_PATH
    output:
        OUTPUT_REP + "02_clustering/clustering_report_full_{fraction}.html"
    script:
        "clustering_report_full.Rmd"
        
#-------------------------------------------------------------------------------
"""
Dummy rule to allow using clusters as wildcard
"""

output = expand(OUTPUT_DAT + "03_sepd/sce_hsc_cluster_{cluster}-sep", cluster = clusters_hsc)
output = output + expand(OUTPUT_DAT + "03_sepd/sce_str_cluster_{cluster}-sep", cluster = clusters_str)
rule separate_sce:
    input: 
        sce_input = expand(OUTPUT_DAT + "02_clst/louvn_clust/sce_{fraction}-02", fraction = fractions)
    output:
        sce_output = output
    script:
        "scripts/03_separate_dummy.R"
        
        
#-------------------------------------------------------------------------------
 
"""
Clustering Annotation
"""

# get marker genes for each cluster 
rule find_markers:
    input: 
        sce_input = rules.louvain_clustering.output
    output:
        markers = OUTPUT_DAT + "04_annc/01_markers/markers_{fraction}"
    params:
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
    script:
        "scripts/04_markers_clusters.R"

# perform preliminary GO analysis for overview (Very basic)
rule go_analysis:
    input: 
        markers = rules.find_markers.output
    output:
        go = OUTPUT_DAT + "04_annc/02_goan/go_{fraction}"
    script:
        "scripts/04_go_clusters.R"

# report on marker gene expression and GO 
rule make_report:
    input: 
        markers = rules.find_markers.output,
        go = rules.go_analysis.output,
        cluster_annotations = "cluster_annotations.txt",
        sce_sep = OUTPUT_DAT + "03_sepd/sce_{fraction}_cluster_{cluster}-sep",
        sce_input = rules.louvain_clustering.output,
        gene_list = "gene_list.txt"
    output:
        OUTPUT_REP + "03_anno_clusters/annotation_{fraction}_cluster_{cluster}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "cluster_anno_report.Rmd"
        
# assign annotation labels informed from report
rule assign_annotation:
    input:
        sce_input = rules.louvain_clustering.output,
        cluster_annotations = "cluster_annotations.txt"
    output:      
        sce_output = OUTPUT_DAT + "04_annc/03_sce/sce_{fraction}-04"
    script:
        "scripts/04_anno_clusters.R"
