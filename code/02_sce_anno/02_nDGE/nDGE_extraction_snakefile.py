#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/02_sce_anno/04_nDGE"

COLORS_REF = config["base"] + config["metadata_paths"]["colors_ref"]
COLORS = config["base"] + config["metadata_paths"]["colors"]

GENES_CLUSTERS = config["base"] + config["metadata_paths"]["gene_list_clusters"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")

clusters_hsc = list(range(1,14))
clusters_hsc = list(map(str,clusters_hsc))

clusters_str = list(range(1,18))
clusters_str = list(map(str,clusters_str))

PADJ_CUTOFF = config["values"]["02_sce_anno"]["ndge_padj_cutoff"]
FC_CUTOFF = config["values"]["02_sce_anno"]["ndge_fc_cutoff"]
tf = "PC_" + str(PADJ_CUTOFF) + "_FC_" + str(FC_CUTOFF)
print(tf)

#-------------------------------------------------------------------------------

# construct paths for all possible outputs/targets, required for rule all
targets = []
for c in clusters_hsc:
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_hsc_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_hsc_cluster_" + c + ".html"] 

for c in clusters_str:
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_str_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_str_cluster_" + c + ".html"] 

for f in fractions:
  targets = targets + [OUTPUT_DAT + "/05_desq/deseq_" + f + "-05"]
  targets = targets + [OUTPUT_DAT + "/06_dsqc/rlog_" + f]
  targets = targets + [OUTPUT_DAT + "/06_dsqc/sva_" + f]
  targets = targets + [OUTPUT_DAT + "/07_tdsq/deseq_" + f + "-07"]
  targets = targets + [OUTPUT_DAT + "/08_nres/" + tf + "/res_" + f + "_species"]
  targets = targets + [OUTPUT_DAT + "/08_nres/" + tf + "/res_" + f + "_cluster"]
  targets = targets + [OUTPUT_DAT + "/08_nres/" + tf + "/res_" + f + "_cluster_dfs"]
  targets = targets + [OUTPUT_DAT + "/08_nres/" + tf + "/res_" + f + "_cluster_shared"]

#print(targets)

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
    cluster="[0-9]+",

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
"""
DESeq pre-processing and QC 
"""

# form pseudobulks/aggregate per cluster per sample and convert to DESeq2 object
# output = list of DESeq2 objects generated per fraction
rule aggregate_convert:
    input:
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04"
    output:
        deseq_output = OUTPUT_DAT + "/05_desq/deseq_{fraction}-05"
    script:
        "scripts/05_aggregate_convert.R"

# produce rlog-transformed values for QC visualisation
# identify hidden sources of variations (SV) for QC of DESeq2 objects
rule qc_deseq:
    input:
        deseq_input = rules.aggregate_convert.output
    output:
        rlog = OUTPUT_DAT + "/06_dsqc/rlog_{fraction}",
        sva = OUTPUT_DAT + "/06_dsqc/sva_{fraction}"
    script:
        "scripts/06_prep_qc_deseq.R"  

# preprocessing/transformation/normalisation of DESeq2 Objects
# SVs are taken into account as batch correction
rule preprocessing_deseq:
    input:
        deseq_input = rules.aggregate_convert.output,
        sva = rules.qc_deseq.output.sva
    output:
        deseq_output = OUTPUT_DAT + "/07_tdsq/deseq_{fraction}-07"
    script:
        "scripts/07_prepro_deseq.R"    

# export results in cluster- or species-specific lists 
rule export_results_ndge:
    input:
        deseq_input = rules.preprocessing_deseq.output
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        species = species,
        clusters_hsc = clusters_hsc,
        clusters_str = clusters_str
    output:
        species_res = OUTPUT_DAT + "/08_nres/" + tf + "/res_{fraction}_species",
        cluster_res = OUTPUT_DAT + "/08_nres/" + tf + "/res_{fraction}_cluster",
        cluster_res_dfs = OUTPUT_DAT + "/08_nres/" + tf + "/res_{fraction}_cluster_dfs",
        cluster_res_list_shared = OUTPUT_DAT + "/08_nres/" + tf + "/res_{fraction}_cluster_shared"
    script:
        "scripts/08_export_ndge_results.R"    
 
#-------------------------------------------------------------------------------
"""   
Make Reports
"""

rule ndge_cluster_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04",
        cluster_res_list_shared = rules.export_results_ndge.output.cluster_res_list_shared,
        cluster_res = rules.export_results_ndge.output.cluster_res,
        gene_list = GENES_CLUSTERS
    output:
        OUTPUT_REP + "/ndge/ndge_report_{fraction}_cluster_{cluster}.html"
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "ndge_cluster_report.Rmd"
  
# this might sometimes fail initially, because this rule doesn't wait long 
# enough for rlog to be saved
rule ndge_bulk_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04",
        rlog = rules.qc_deseq.output.rlog
    output:
        OUTPUT_REP + "/bulk/bulk_quality_report_{fraction}_cluster_{cluster}.html"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "ndge_bulk_quality_report.Rmd"
