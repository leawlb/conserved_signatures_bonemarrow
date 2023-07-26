#!/bin/python 

#!/bin/python 

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno/"
OUTPUT_REP = OUTPUT_BASE + "/reports/02_sce_anno/05_subclustering/"

# clusters for subclustering
# hsc clusters 2 and 4 will be merged later
clusters_hsc = ["2", "6"]
clusters_str = ["4"]

subclusters_hsc = list(range(1,13))
subclusters_hsc = list(map(str,subclusters_hsc))

subclusters_str = list(range(1,9))
subclusters_str = list(map(str,subclusters_str))

#-------------------------------------------------------------------------------

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")

# collect all paths for all possible outputs/targets, required for rule all
targets = []
for c in clusters_hsc:
  targets = targets + [OUTPUT_DAT + "09_sbcl/sce_hsc_cluster_" + c + "-09"]

for c in clusters_str:
  targets = targets + [OUTPUT_DAT + "09_sbcl/sce_str_cluster_" + c + "-09"]

for f in fractions:
  targets = targets + [OUTPUT_DAT + "10_anns/sce_" + f + "-10"]
  targets = targets + [OUTPUT_REP + "subclustering_report_" + f + ".html"]
  
for s in subclusters_hsc:
  targets = targets + [OUTPUT_DAT + "11_sepc/sce_hsc_subcluster_" + s + "-sep"]
  #targets = targets + [OUTPUT_REP + "03_anno_clusters/annotation_hsc_cluster_" + c + ".html" ] 

for s in subclusters_str:
  targets = targets + [OUTPUT_DAT + "11_sepc/sce_str_subcluster_" + s + "-sep"]
  #targets = targets + [OUTPUT_REP + "03_anno_clusters/annotation_str_cluster_" + c + ".html"] 

# local execution of non-demanding rules
localrules: all  
 
rule all: # must contain all possible output paths from all rules
    input:
        targets
        
#-------------------------------------------------------------------------------
# subclustering
rule subclustering:
    input:
        sce_sep = OUTPUT_DAT + "03_sepd/sce_{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "04_annc/03_sce/sce_{fraction}-04",
        subclustering_genes = "subclustering_genes.txt",
        subclustering_annotation = "subclustering_annotation.txt",
        genes_shared_list = OUTPUT_DAT + "08_dres/PC_0.05_FC_1.5/res_{fraction}_cluster_shared"
    output:
        sce_output = OUTPUT_DAT + "09_sbcl/sce_{fraction}_cluster_{cluster}-09"
    script:
        "scripts/09_subclustering.R" 
  
# add subclusters to fraction object and annotate
sce_subcl_input = expand(rules.subclustering.output, cluster = clusters_hsc, fraction = ["hsc"]) 
sce_subcl_input = sce_subcl_input + expand(rules.subclustering.output, cluster = clusters_str, fraction = ["str"])
print(sce_subcl_input)
rule add_subclusters:
    input:
        sce_input = OUTPUT_DAT + "04_annc/03_sce/sce_{fraction}-04",
        sce_subcl = sce_subcl_input,
        final_annotation = "final_annotation.txt"
    output:
        sce_output = OUTPUT_DAT + "10_anns/sce_{fraction}-10"
    script:
        "scripts/10_add_subclusters_anno.R" 

# report
rule make_subclustering_report:
    input:
        sce_input = rules.add_subclusters.output,
        gene_list_subclustering = "subclustering_genes.txt"
    params:
        color_table = "colors.txt"
    output:
        OUTPUT_REP + "subclustering_report_{fraction}.html"
    script:
        "subclustering_report.Rmd" 
  
"""   
Dummy rule to allow using subclusters as wildcard
"""

output = expand(OUTPUT_DAT + "11_sepc/sce_hsc_subcluster_{subcluster}-sep", subcluster = subclusters_hsc)
output = output + expand(OUTPUT_DAT + "11_sepc/sce_str_subcluster_{subcluster}-sep", subcluster = subclusters_str)
print(output)
rule separate_sce:
    input: 
        sce_input = expand(OUTPUT_DAT + "10_anns/sce_{fraction}-10", fraction = fractions)
    output:
        sce_output = output
    script:
        "scripts/11_separate_dummy.R"
 
"""     
# repeat marker gene and GO analysis on the subclusters now
# get marker genes for each cluster 
rule find_markers:
    input: 
        sce_input = rules.separate_sce.output
    output:
        markers = OUTPUT_DAT + "11_anqc/01_markers/markers_{fraction}"
    params:
        nr_hvgs = config["values"]["preprocessing"]["nr_hvgs"],
    script:
        "scripts/11_markers_sublusters.R"

# perform preliminary GO analysis for overview (Very basic)
rule go_analysis:
    input: 
        markers = rules.find_markers.output
    output:
        go = OUTPUT_DAT + "11_anqc/02_goan/go_{fraction}"
    script:
        "scripts/11_go_subclusters.R"

# report on marker gene expression and GO 
rule make_report_subcl_markers:
    input: 
        markers = rules.find_markers.output,
        go = rules.go_analysis.output,
        cluster_annotations = "cluster_annotations.txt",
        sce_sep = OUTPUT_DAT + "11_sepc/sce_{fraction}_subcluster_{subcluster}-sep",
        sce_input = rules.add_subclusters.output,
        gene_list = "gene_list.txt"
    output:
        OUTPUT_REP + "subclustering_qc/annotation_{fraction}_cluster_{cluster}.html"
    params:
        color_tables = TABLES_PATH
    script:
        "subclustering_qc_report.Rmd"

"""
