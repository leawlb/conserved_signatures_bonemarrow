#!/bin/python 

"""
Batch correction, clustering and cluster annotation using louvain clustering. 
"""

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/02_sce_anno"

ANNO_CLUSTERS = config["base"] + config["metadata_paths"]["annotation_clusters"]
GENES_CLUSTERS = config["base"] + config["metadata_paths"]["gene_list_clusters"]

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
COLORS = config["base"] + config["metadata_paths"]["colors"]

VALUES =  config["values"]["02_sce_anno"]
BATCH_USE = VALUES["batch_use"] # which colData to use as batch

RUN_ANNO_REPORTS = config["run_anno_reports"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")
ages = get_list(metadata = METADATA, column = "Age_ID")

# cannot take from metadata since cluster number is dependent on fraction
clusters_hsc = list(range(1,14))
clusters_hsc = list(map(str, clusters_hsc))

clusters_str = list(range(1,21))
clusters_str = list(map(str, clusters_str))

print(clusters_hsc)
print(clusters_str)

#-------------------------------------------------------------------------------

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_mnnc/sce_" + f + "_" + BATCH_USE + "-01"]
  targets = targets + [OUTPUT_DAT + "/02_clst/sce_" + f + "-02"]
  targets = targets + [OUTPUT_REP + "/01_batch_correction/batch_correction_report_" + f + "_" + BATCH_USE + ".html"]
  targets = targets + [OUTPUT_REP + "/02_clustering/clustering_report_" + f +".html"]
  targets = targets + [OUTPUT_DAT + "/04_annc/01_markers/markers_" + f]
  targets = targets + [OUTPUT_DAT + "/04_annc/02_goan/go_" + f]
  targets = targets + [OUTPUT_DAT + "/04_annc/03_sce/sce_" + f + "-04"]
# 
for c in clusters_hsc:
  targets = targets + [OUTPUT_DAT + "/03_sepd/hsc_cluster_" + c + "-sep"]
  if RUN_ANNO_REPORTS:
    targets = targets + [OUTPUT_REP + "/03_anno_clusters/annotation_hsc_cluster_" + c + ".html" ]
 
for c in clusters_str:
   targets = targets + [OUTPUT_DAT + "/03_sepd/str_cluster_" + c + "-sep"]
   if RUN_ANNO_REPORTS:
    targets = targets + [OUTPUT_REP + "/03_anno_clusters/annotation_str_cluster_" + c + ".html"]

for s in species:
  for f in fractions:
    targets = targets + [OUTPUT_DAT + "/01_mnnc/comparison/sce_" + s + "_" + f + "_" + BATCH_USE + "-01"]
    targets = targets + [OUTPUT_DAT + "/02_clst/comparison/sce_" + s + "_" + f + "-02"]
    targets = targets + [OUTPUT_REP + "/02_clustering/comparison/clustering_report_" + s + "_" + f + ".html"]
  
#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------
"""
Info on batch correction methods from:
Luecken et al. "Benchmarking atlas-level data integration in single-cell 
genomics", Nat Met 2022
Tran, Ang, Chevrier, Zhang et al. "A benchmark of batch effect correction 
methods for single-cell RNA sequencing data", Genome Biology 2020

Batch correction using MNNcorrect (fastMNN):

"""
rule run_mnncorrect:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/fractions/sce_{fraction}-08"
    output:
        sce_output = OUTPUT_DAT + "/01_mnnc/sce_{fraction}_" + BATCH_USE + "-01"
    params:
        batch_use = BATCH_USE,
        nr_hvgs_BC = VALUES["nr_hvgs_batch_correction"], 
        nr_hvgs = config["values"]["nr_hvgs"], 
        seeds_umap = VALUES["seeds_umap"]
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
        OUTPUT_REP + "/01_batch_correction/batch_correction_report_{fraction}_" + BATCH_USE + ".html"
    params:
        batch_use = BATCH_USE,
        colors_ref_path = COLORS_REF,
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
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
        sce_output = OUTPUT_DAT + "/02_clst/sce_{fraction}-02"
    params:
        k_graph_list = VALUES["k_graph_list"],
        resolution_louvain_list = VALUES["resolution_louvain_list"]
    script:
        "scripts/02_louvain_clustering.R"
  
"""
Make clustering reports 
"""      
        
rule make_clustering_report:
    input:
        sce_l = rules.louvain_clustering.output
    params:
        colors_ref_path = COLORS_REF,
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    output:
        OUTPUT_REP + "/02_clustering/clustering_report_{fraction}.html"
    script:
        "clustering_report.Rmd"
        
#-------------------------------------------------------------------------------
"""
This only saves the cluster number, not a SCE object, but it allows using
clusters as a wildcard
"""

output = expand(OUTPUT_DAT + "/03_sepd/hsc_cluster_{cluster}-sep", cluster = clusters_hsc)
output = output + expand(OUTPUT_DAT + "/03_sepd/str_cluster_{cluster}-sep", cluster = clusters_str)
rule separate_sce:
    input: 
        sce_input = expand(OUTPUT_DAT + "/02_clst/sce_{fraction}-02", fraction = fractions)
    output:
        output = output
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
        markers = OUTPUT_DAT + "/04_annc/01_markers/markers_{fraction}"
    params:
        nr_hvgs = config["values"]["nr_hvgs"]
    script:
        "scripts/04_markers_clusters.R"

# perform preliminary GO analysis for overview (very basic)
rule go_analysis:
    input: 
        markers = rules.find_markers.output
    output:
        go = OUTPUT_DAT + "/04_annc/02_goan/go_{fraction}"
    script:
        "scripts/04_go_clusters.R"

if RUN_ANNO_REPORTS:
  # report on marker gene expression and GO 
  rule make_anno_report:
      input: 
          markers = rules.find_markers.output,
          go = rules.go_analysis.output,
          sep = OUTPUT_DAT + "/03_sepd/{fraction}_cluster_{cluster}-sep",
          sce_input = rules.louvain_clustering.output,
          gene_list = GENES_CLUSTERS
      output:
          OUTPUT_REP + "/03_anno_clusters/annotation_{fraction}_cluster_{cluster}.html"
      params:
          colors_ref_path = COLORS_REF,
          colors_path = COLORS,
          functions = "../../source/sce_functions.R",
          plotting = "../../source/plotting.R",
          colors = "../../source/colors.R"
      script:
          "cluster_anno_report.Rmd"
        
# assign annotation labels informed from report
rule assign_annotation:
    input:
        sce_input = rules.louvain_clustering.output,
        anno_clusters = ANNO_CLUSTERS
    output:      
        sce_output = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04"
    script:
        "scripts/04_anno_clusters.R"
        
#-------------------------------------------------------------------------------
"""
Cluster species separately and compare to original clustering
Re-use the same scripts
"""

# batch correction
rule run_mnncorrect_species:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/01_sce_prep/08_mrge/species/sce_{species}_{fraction}-08"
    output:
        sce_output = OUTPUT_DAT + "/01_mnnc/comparison/sce_{species}_{fraction}_" + BATCH_USE + "-01"
    params:
        batch_use = BATCH_USE,
        nr_hvgs_BC = VALUES["nr_hvgs_batch_correction"], 
        nr_hvgs = config["values"]["nr_hvgs"], 
        seeds_umap = VALUES["seeds_umap"]
    script:
        "scripts/01_mnncorrect.R"

# clustering
rule louvain_clustering_species:
    input: 
        sce_input = rules.run_mnncorrect_species.output
    output:
        sce_output = OUTPUT_DAT + "/02_clst/comparison/sce_{species}_{fraction}-02"
    params:
        k_graph_list = VALUES["k_graph_list"],
        resolution_louvain_list = VALUES["resolution_louvain_list"]
    script:
        "scripts/02_louvain_clustering.R"

"""
Make clustering comparison reports 
"""    
        
rule make_clustering_comparison_report:
    input:
        sce_l = rules.assign_annotation.output,
        sce_l_species = rules.louvain_clustering_species.output,
        sce_l_celltypes = rules.assign_annotation.output
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    conda:
        "../../envs/vis_ggalluvial.yml"
    output:
        OUTPUT_REP + "/02_clustering/comparison/clustering_report_{species}_{fraction}.html"
    script:
        "clustering_comparison_report.Rmd"
