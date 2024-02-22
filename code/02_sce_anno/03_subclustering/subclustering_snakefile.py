#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/02_sce_anno/05_subclustering"

COLORS = config["base"] + config["metadata_paths"]["colors"]
COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]

VALUES =  config["values"]["02_sce_anno"]

GENE_LIST_CLUSTERS = config["base"] + config["metadata_paths"]["gene_list_clusters"]
GENE_LIST_SUBCLUSTERING = config["base"] + config["metadata_paths"]["gene_list_subclustering"]
GENE_LIST_DOTPLOT = config["base"] + config["metadata_paths"]["gene_list_dotplot"]
ANNO_SUBCLUSTERS = config["base"] + config["metadata_paths"]["annotation_subclusters"]
ANNO_FINAL = config["base"] + config["metadata_paths"]["annotation_final"]

#-------------------------------------------------------------------------------

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")

# clusters for subclustering
# hsc clusters 2 and 4 will be merged later
clusters_hsc = ["2", "7"]
clusters_str = ["4"]

# subclusters for dummy separation
subclusters_hsc = list(range(1,13))
subclusters_hsc = list(map(str,subclusters_hsc))

subclusters_str = list(range(1,9))
subclusters_str = list(map(str,subclusters_str))

#-------------------------------------------------------------------------------

targets = []
for c in clusters_hsc:
  targets = targets + [OUTPUT_DAT + "/09_sbcl/sce_hsc_cluster_" + c + "-09"]
  targets = targets + [OUTPUT_REP + "/mclust/mclust_report_hsc_cluster_" + c + ".html"]

for c in clusters_str:
  targets = targets + [OUTPUT_DAT + "/09_sbcl/sce_str_cluster_" + c + "-09"]
  targets = targets + [OUTPUT_REP + "/mclust/mclust_report_str_cluster_" + c + ".html"]

for f in fractions:
  targets = targets + [OUTPUT_DAT + "/10_anns/sce_" + f + "-10"]
  targets = targets + [OUTPUT_REP + "/results/results_report_" + f + ".html"]
  targets = targets + [OUTPUT_DAT + "/12_anqc/01_markers/markers_" + f]
  targets = targets + [OUTPUT_DAT + "/12_anqc/02_go/go_" + f]

for s in subclusters_hsc:
  targets = targets + [OUTPUT_DAT + "/11_sepc/sce_hsc_subcluster_" + s + "-sep"]
  targets = targets + [OUTPUT_REP + "/markergenes/markergenes_report_hsc_subcluster_" + s + ".html"]
  
for s in subclusters_str:
  targets = targets + [OUTPUT_DAT + "/11_sepc/sce_str_subcluster_" + s + "-sep"]
  targets = targets + [OUTPUT_REP + "/markergenes/markergenes_report_str_subcluster_" + s + ".html"]

#-------------------------------------------------------------------------------

localrules: all  
 
rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
# subclustering
rule subclustering:
    input:
        sce_sep = OUTPUT_DAT + "/03_sepd/sce_{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04",
        gene_list_subcl = GENE_LIST_SUBCLUSTERING,
        anno_subcl = ANNO_SUBCLUSTERS,
        genes_list_shared = OUTPUT_DAT + "/08_nres/PC_0.05_FC_1.5/res_{fraction}_cluster_shared"
    output:
        sce_output = OUTPUT_DAT + "/09_sbcl/sce_{fraction}_cluster_{cluster}-09"
    script:
        "scripts/09_subclustering.R" 
        
# report on subclustering genes, mclust, subclusters per cluster
rule subclustering_mclust_report:
    input:
        sce_input = rules.subclustering.output, #
        gene_list_subclustering = GENE_LIST_SUBCLUSTERING, #
        genes_list_shared = OUTPUT_DAT + "/08_nres/PC_0.05_FC_1.5/res_{fraction}_cluster_shared" #
    output:
        OUTPUT_REP + "/mclust/mclust_report_{fraction}_cluster_{cluster}.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "subclustering_mclust_report.Rmd" 
  
# add subclusters to fraction object and annotate
sce_subcl_input = expand(rules.subclustering.output, cluster = clusters_hsc, fraction = ["hsc"]) 
sce_subcl_input = sce_subcl_input + expand(rules.subclustering.output, cluster = clusters_str, fraction = ["str"])
print(sce_subcl_input)
rule add_subclusters:
    input:
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04",
        sce_subcl = sce_subcl_input,
        anno_final = ANNO_FINAL
    output:
        sce_output = OUTPUT_DAT + "/10_anns/sce_{fraction}-10"
    script:
        "scripts/10_add_subclusters_anno.R" 

# report results of clustering for entire fraction
rule make_subclustering_results_report:
    input:
        sce_input = rules.add_subclusters.output,
        gene_list_subclustering = GENE_LIST_SUBCLUSTERING,
        gene_list_dtplt = GENE_LIST_DOTPLOT,
        genes_list_shared = OUTPUT_DAT + "/08_nres/PC_0.05_FC_1.5/res_{fraction}_cluster_shared"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    output:
        OUTPUT_REP + "/results/results_report_{fraction}.html"
    script:
        "subclustering_results_report.Rmd" 
  
# Dummy rule to allow using subclusters as wildcard
output = expand(OUTPUT_DAT + "/11_sepc/sce_hsc_subcluster_{subcluster}-sep", subcluster = subclusters_hsc)
output = output + expand(OUTPUT_DAT + "/11_sepc/sce_str_subcluster_{subcluster}-sep", subcluster = subclusters_str)
print(output)
rule separate_sce:
    input: 
        sce_input = expand(OUTPUT_DAT + "/10_anns/sce_{fraction}-10", fraction = fractions)
    output:
        sce_output = output
    script:
        "scripts/11_separate_dummy.R"
 
# repeat marker gene and GO analysis on the subclusters now
# get marker genes for each cluster 
rule find_markers:
    input: 
        sce_input = rules.add_subclusters.output
    output:
        markers = OUTPUT_DAT + "/12_anqc/01_markers/markers_{fraction}"
    params:
        nr_hvgs = config["values"]["nr_hvgs"],
    script:
        "scripts/12_markers_subclusters.R"

# perform preliminary GO analysis for overview (Very basic)
rule go_analysis:
    input: 
        markers = rules.find_markers.output
    output:
        go = OUTPUT_DAT + "/12_anqc/02_go/go_{fraction}"
    script:
        "scripts/12_go_subclusters.R"

# report on marker gene expression and GO 
rule make_subclustering_markers_report:
    input: 
        markers = rules.find_markers.output,
        go = rules.go_analysis.output,
        sce_sep = OUTPUT_DAT + "/11_sepc/sce_{fraction}_subcluster_{subcluster}-sep",
        sce_input = rules.add_subclusters.output,
        gene_list_subcl = GENE_LIST_CLUSTERS
    output:
        OUTPUT_REP + "/markergenes/markergenes_report_{fraction}_subcluster_{subcluster}.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "subclustering_markergenes_report.Rmd"
