#!/bin/python 

import pandas as pd

OUTPUT_BASE = config["paths"]["output_dir"]
METADATA = pd.read_csv(config["metadata"]["table"])

# specific data and report output paths
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/02_sce_anno/"
OUTPUT_REP = OUTPUT_BASE + "/reports/02_sce_anno/04_nDGE/"

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
species = get_list(metadata = METADATA, column = "Species_ID")

padj_cutoff = config["values"]["annotation"]["ndge_padj_cutoff"]
fc_cutoff = config["values"]["annotation"]["ndge_fc_cutoff"]
tf = "PC_" + str(padj_cutoff) + "_FC_" + str(fc_cutoff)
print(tf)



# construct paths for all possible outputs/targets, required for rule all
targets = []
for c in clusters_hsc:
  targets = targets + [OUTPUT_REP + "ndge_report_hsc_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "bulk_quality_report_hsc_cluster_" + c + ".html"] 

for c in clusters_str:
  targets = targets + [OUTPUT_REP + "ndge_report_str_cluster_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "bulk_quality_report_str_cluster_" + c + ".html"] 

for f in fractions:
  targets = targets + [OUTPUT_DAT + "05_desq/deseq_" + f + "-05"]
  targets = targets + [OUTPUT_DAT + "06_dsqc/rlog_" + f]
  targets = targets + [OUTPUT_DAT + "06_dsqc/sva_" + f]
  targets = targets + [OUTPUT_DAT + "07_tdsq/deseq_" + f + "-07"]
  targets = targets + [OUTPUT_DAT + "08_dres/" + tf + "/res_" + f + "_species"]
  targets = targets + [OUTPUT_DAT + "08_dres/" + tf + "/res_" + f + "_cluster"]
  targets = targets + [OUTPUT_DAT + "08_dres/" + tf + "/res_" + f + "_cluster_dfs"]
  targets = targets + [OUTPUT_DAT + "08_dres/" + tf + "/res_" + f + "_cluster_shared"]

#print(targets)

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
    cluster="[0-9]+",

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
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
        sce_input = OUTPUT_DAT + "04_annc/03_sce/sce_{fraction}-04"
    output:
        deseq_output = OUTPUT_DAT + "05_desq/deseq_{fraction}-05"
    script:
        "scripts/05_aggregate_convert.R"

# produce rlog-transformed values for QC visualisation
# identify hidden sources of variations (SV) for QC of DESeq2 objects
rule qc_deseq:
    input:
        deseq_input = rules.aggregate_convert.output
    output:
        rlog = OUTPUT_DAT + "06_dsqc/rlog_{fraction}",
        sva = OUTPUT_DAT + "06_dsqc/sva_{fraction}"
    script:
        "scripts/06_prep_qc_deseq.R"  

# preprocessing/transformation/normalisation of DESeq2 Objects
# SVs are taken into account as batch correction
rule preprocessing_deseq:
    input:
        deseq_input = rules.aggregate_convert.output,
        sva = rules.qc_deseq.output.sva
    output:
        deseq_output = OUTPUT_DAT + "07_tdsq/deseq_{fraction}-07"
    script:
        "scripts/07_prepro_deseq.R"    

# export results in cluster- or species-specific lists 
rule export_results_ndge:
    input:
        deseq_input = rules.preprocessing_deseq.output
    params:
        padj_cutoff = padj_cutoff,
        fc_cutoff = fc_cutoff,
        species = species,
        clusters_hsc = clusters_hsc,
        clusters_str = clusters_str
    output:
        species_res = OUTPUT_DAT + "08_dres/" + tf + "/res_{fraction}_species",
        cluster_res = OUTPUT_DAT + "08_dres/" + tf + "/res_{fraction}_cluster",
        cluster_res_dfs = OUTPUT_DAT + "08_dres/" + tf + "/res_{fraction}_cluster_dfs",
        cluster_res_list_shared = OUTPUT_DAT + "08_dres/" + tf + "/res_{fraction}_cluster_shared"
    script:
        "scripts/08_export_ndge_results.R"    
 

"""   

Make Reports

"""

rule cluster_ndge_report:
    input: 
        sce_sep = OUTPUT_DAT + "03_sepd/sce_{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "04_annc/03_sce/sce_{fraction}-04",
        cluster_res_list_shared = rules.export_results_ndge.output.cluster_res_list_shared,
        cluster_res = rules.export_results_ndge.output.cluster_res,
        gene_list = "../01_clustering/gene_list.txt"
    output:
        OUTPUT_REP + "ndge_report_{fraction}_cluster_{cluster}.html"
    params:
        padj_cutoff = padj_cutoff,
        fc_cutoff = fc_cutoff
    script:
        "make_cluster_ndge_report.Rmd"
  



rule cluster_bulk_report:
    input: 
        sce_sep = OUTPUT_DAT + "03_sepd/sce_{fraction}_cluster_{cluster}-sep",
        sce_input = OUTPUT_DAT + "04_annc/03_sce/sce_{fraction}-04",
        rlog = rules.qc_deseq.output.rlog
        #sva_11 = rules.qc_dge.output.sva_11,
        #tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        #tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11
    output:
        OUTPUT_REP + "bulk_quality_report_{fraction}_cluster_{cluster}.html"
    script:
        "make_bulk_quality_report.Rmd"

"""

# report on shared genes for pval and cutoff selection
rule make_report_shared:
    input: 
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/dsq_objects/tdsq_after_{fraction}-11", # after batch correction
        cluster_res_df = rules.export_results_ndge.output.cluster_res_df,
        gene_list = config["metadata"]["gene_list_hsc"],
    params:
        padj_cutoff = config["values"]["DGE"]["flayer_integrated_padj_cutoff"],
        fc_cutoff = config["values"]["DGE"]["flayer_integrated_fc_cutoff"]
    output:
        OUTPUT_BASE_PATH + "/reports/009_DGE_analysis_flayer/{fraction}/shared_nDGE/" + tf + "/shared_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_shared.Rmd"







rule export_results_visualisation_before:
    input:
        deseq_input = rules.preprocessing_deseq.output
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/res_objects/DGE/res_before_{fraction}-11"
    script:
        "scripts/12_export_results_vis_flayer_integrated.R"    





rule make
     
     
   
# get results from all pairwise comparisons
rule export_results_visualisation_before:
    input:
        deseq_input = rules.preprocessing_deseq.output
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/res_objects/DGE/res_before_{fraction}-11"
    script:
        "scripts/12_export_results_vis_flayer_integrated.R"    



# progress report on pseudo-bulk quality, SVs and effect of BC with SVs
rule make_report_flayer_qc_bulks:
    input: 
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        rld_11 = rules.qc_dge.output.rld_11,
        sva_11 = rules.qc_dge.output.sva_11,
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11
    output:
        OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/{fraction}/quality_bulks/quality_bulks_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_qc_dge.Rmd"



# get results from pairwise comparisons for visualisation
rule export_results_visualisation_after:
    input:
        tdsq = rules.preprocessing_dge.output.tdsq_after_11
    output:
        res = OUTPUT_BASE_PATH + "/sce_objects/12_dge_prepro/res_objects/DGE/res_after_{fraction}-11"
    script:
        "scripts/12_export_results_vis_flayer_integrated.R"    

# progress report on DGE quality of DGE corrected for batch, age, SVs
rule make_report_flayer_pairwise_dge:
    input: 
        sce_11 = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_anno_{fraction}-11",
        sce_sep = OUTPUT_BASE_PATH + "/sce_objects/11_annotation/sce_objects/sce_{fraction}_cluster_{cluster}-sep",
        tdsq_before_11 =rules.preprocessing_dge.output.tdsq_before_11,
        tdsq_after_11 = rules.preprocessing_dge.output.tdsq_after_11,
        res_before_11 =rules.export_results_visualisation_before.output,
        res_after_11 = rules.export_results_visualisation_after.output
    output:
        OUTPUT_BASE_PATH + "/reports/008_DGE_preprocessing/{fraction}/quality_DGE/quality_dge_{fraction}_cluster_{cluster}.html"
    script:
        "report_flayer_pairwise_dge.Rmd"
   
"""
