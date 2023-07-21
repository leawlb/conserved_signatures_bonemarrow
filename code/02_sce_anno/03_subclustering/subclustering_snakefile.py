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
  
# local execution of non-demanding rules
localrules: all  
 
rule all: # must contain all possible output paths from all rules
    input:
        targets
        
#-------------------------------------------------------------------------------
# subclustering
"""
!!!! FOR THIS STEP, ANOTHER ENVIRONMENT IS REQUIRED:
use snakemake_isbm_mclust
USE THAT ENVIRONMENT FOR THIS STEP ONLY!!
!!!!
"""
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
sce_subcl_input = expand(rules.subclustering.output, cluster = clusters_hsc, fraction = ["hsc"]) + expand(rules.subclustering.output, cluster = clusters_str, fraction = ["str"])
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

"""
# report
rule make_subclustering_report:
    input:
        sce_15 = rules.anno_threelayers.output.sce_15_anno,
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        gene_list_subclustering = config["metadata"]["path"] + "/gene_list_{fraction}_subclustering.csv"
    output:
        OUTPUT_BASE_PATH + "/reports/011_slayer_annotation/{fraction}/annotated/anno_report_{fraction}_{cluster}.html"
    script:
        "report_anno.Rmd" 
        
# summary

rule make_anno_summary:
    input:
        sce_15 = rules.anno_threelayers.output,
        gene_list_subclustering = config["metadata"]["path"] + "/gene_list_{fraction}_subclustering.csv"
    output:
        OUTPUT_BASE_PATH + "/reports/011_slayer_annotation/{fraction}/annotated/anno_summary_{fraction}.html"
    script:
        "summary_anno.Rmd" 
"""
