#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/02_sce_anno/04_nDGE"

COLORS = config["base"] + config["metadata_paths"]["colors"]

GENES_CLUSTERS = config["base"] + config["metadata_paths"]["gene_list_clusters"]
SV_PATH = config["base"] + config["metadata_paths"]["sources_variation"]["annotation_species"]

METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")

# only need to do nDGE analysis here for clusters to be subclustered
# for all other clusters, nDGE is performed in 03_sce_analysis
clusters_hsc = ["2", "4", "6"]
clusters_str = ["6"]

PADJ_CUTOFF = config["values"]["02_sce_anno"]["ndge_padj_cutoff"]
FC_CUTOFF = config["values"]["02_sce_anno"]["ndge_fc_cutoff"]
tf = "PC_" + str(PADJ_CUTOFF) + "_FC_" + str(FC_CUTOFF)
print(tf)

#-------------------------------------------------------------------------------

targets = []
for c in clusters_hsc:
  #targets = targets + [OUTPUT_REP + "/ndge/ndge_report_cluster_" + c + "_hsc.html"] 
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_sv_" + c + "_hsc.html"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_cluster_" + c + "_hsc.html"] 

for c in clusters_str:
  #targets = targets + [OUTPUT_REP + "/ndge/ndge_report_cluster_" + c + "_str.html"] 
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_sv_" + c + "_str.html"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_cluster_" + c + "_str.html"] 

for f in fractions:
  targets = targets + [OUTPUT_DAT + "/05_desq/deseq_" + f + "-05"]
  targets = targets + [OUTPUT_DAT + "/06_dsqc/rlog_" + f]
  targets = targets + [OUTPUT_DAT + "/06_dsqc/sva_" + f]
  targets = targets + [OUTPUT_DAT + "/07_tdsq/deseq_" + f + "-07"]
  targets = targets + [OUTPUT_DAT + "/08_nres/" + tf + "/res_" + f + "_cluster"]
  targets = targets + [OUTPUT_DAT + "/08_nres/" + tf + "/res_" + f + "_cluster_dfs"]
  targets = targets + [OUTPUT_DAT + "/08_nres/" + tf + "/res_" + f + "_cluster_shared"]

print(targets)

#-------------------------------------------------------------------------------

wildcard_constraints:
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
    resources:
        mem_mb=25000,
        queue = "medium"
    threads: 4
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
    resources:
        mem_mb=5000,
        queue = "medium"
    threads: 4
    script:
        "scripts/06_prep_qc_deseq.R"  

# preprocessing/transformation/normalisation of DESeq2 Objects
# SVs are taken into account as batch correction after evaluation in 
# cluster_sv report
rule preprocessing_deseq:
    input:
        deseq_input = rules.aggregate_convert.output,
        sva = rules.qc_deseq.output.sva,
        sv_path = SV_PATH
    resources:
        mem_mb=5000,
        queue = "medium"
    threads: 4
    output:
        deseq_output = OUTPUT_DAT + "/07_tdsq/deseq_{fraction}-07"
    script:
        "scripts/07_prepro_deseq.R"    

# export results in cluster-specific lists 
rule export_results_ndge:
    input:
        deseq_input = rules.preprocessing_deseq.output
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        species = species,
        clusters_hsc = clusters_hsc,
        clusters_str = clusters_str
    resources:
        mem_mb=5000,
        queue = "medium"
    threads: 4
    output:
        # lists sorted by cluster 
        cluster_res = OUTPUT_DAT + "/08_nres/" + tf + "/res_{fraction}_cluster",
        cluster_resdf_list = OUTPUT_DAT + "/08_nres/" + tf + "/res_{fraction}_cluster_dfs",
        cluster_sharedgenes_list = OUTPUT_DAT + "/08_nres/" + tf + "/res_{fraction}_cluster_shared"
    script:
        "scripts/08_export_ndge_results.R"    
 
#-------------------------------------------------------------------------------
"""   
Make Reports
"""
# check cluster bulk quality
rule ndge_bulk_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04",
        rlog = rules.qc_deseq.output.rlog
    output:
        OUTPUT_REP + "/bulk/bulk_quality_report_cluster_{cluster}_{fraction}.html"
    resources:
        mem_mb=25000,
        queue = "medium"
    threads: 4
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "ndge_bulk_quality_report.Rmd"
  
# check hidden sources of variations and decide which ones to add to DESeq2
# design
rule ndge_sv_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04",
        dsq_list = rules.aggregate_convert.output,
        sva_list = rules.qc_deseq.output.sva
    output:
        OUTPUT_REP + "/ndge/ndge_report_sv_{cluster}_{fraction}.html"
    resources:
        mem_mb=25000,
        queue = "medium"
    threads: 4
    params:
        plotting = "../../source/plotting.R"
    script:
        "ndge_sv_report.Rmd"


# check nDGEs afterwards   
# evaluate which nDGEs (shared by three of four species) would be useful 
# for subclustering
rule ndge_cluster_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "/04_annc/03_sce/sce_{fraction}-04",
        cluster_res = rules.export_results_ndge.output.cluster_res,
        cluster_sharedgenes_list = rules.export_results_ndge.output.cluster_sharedgenes_list,
        gene_list = GENES_CLUSTERS
    output:
        OUTPUT_REP + "/ndge/ndge_report_cluster_{cluster}_{fraction}.html"
    resources:
        mem_mb=15000,
        queue = "medium"
    threads: 4
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "ndge_cluster_report.Rmd"

