#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_conserved_EMF"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_conserved_EMF"

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
print(METADATA)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")
species = get_list(metadata = METADATA, column = "Species_ID")

CELLTYPES = pd.read_csv(config["base"] + config["metadata_paths"]["annotation_final"], sep = ";")

# make cell types list from metadata
values = CELLTYPES.iloc[:,1]
values = values.drop_duplicates()
values = values.squeeze()
values = values.tolist()
celltypes = values
converter = lambda x: x.replace(' ', '_')
celltypes = list(map(converter, celltypes))
converter = lambda x: x.replace('/', '_')
celltypes = list(map(converter, celltypes))
converter = lambda x: x.replace('.', '')
celltypes = list(map(converter, celltypes))
celltypes_hsc = celltypes[0:12]
celltypes_str = celltypes[12:20]
print(celltypes_hsc)
print(celltypes_str)

PADJ_CUTOFF = config["values"]["02_sce_anno"]["ndge_padj_cutoff"]
FC_CUTOFF = config["values"]["02_sce_anno"]["ndge_fc_cutoff"]
tf = "PC_" + str(PADJ_CUTOFF) + "_FC_" + str(FC_CUTOFF)
print(tf)

#-------------------------------------------------------------------------------

# construct paths for all possible outputs/targets, required for rule all
targets = []

for c in celltypes_hsc:  
  targets = targets + [OUTPUT_REP + "/celltypes/emf_celltype_report_hsc_" + c + ".html"] 

for c in celltypes_str:
  targets = targets + [OUTPUT_REP + "/celltypes/emf_celltype_report_str_" + c + ".html"] 

for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_emfs/emf_list_" + f]
  targets = targets + [OUTPUT_DAT + "/02_clst/sce_" + f]
  targets = targets + [OUTPUT_REP + "/clustering_ctrl_report_" + f + ".html"] 
  targets = targets + [OUTPUT_DAT + "/03_cref/sce_baccin_" + f]
  targets = targets + [OUTPUT_DAT + "/03_cref/sce_ref2_" + f]
#
targets = targets + [OUTPUT_REP + "/emf_summary.html"]

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
    fraction="[a-z]+"

localrules: all  

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------

"""   
# Get overlap with Veronica's genes with conserved marker gene function
# cons_markers were provided by Veronica Busa

# These genes are "genes with conserved Expression level and Marker Function"
# = EMF
"""

rule export_emf_genes:
    input:
        celltype_res_list_shared = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/05_nres/PC_0.05_FC_1.5/res_{fraction}_celltype_shared",
        marker_cons = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_{fraction}.RData"
    output:
        cons_EMF_list = OUTPUT_DAT + "/01_emfs/emf_list_{fraction}"
    script:
        "scripts/01_export_emfs.R"    

rule emf_celltype_report:
    input: 
        sep = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies" + "/06_sepd/{fraction}_{celltype}-sep",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        cons_EMF_list = rules.export_emf_genes.output,
        celltype_res_list_shared = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/05_nres/PC_0.05_FC_1.5/res_{fraction}_celltype_shared",
    output:
        OUTPUT_REP + "/celltypes/emf_celltype_report_{fraction}_{celltype}.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "emf_celltype_report.Rmd"

    
rule core_summary:
    input: 
        emf_list_hsc = OUTPUT_DAT + "/01_emfs/emf_list_hsc",
        emf_list_str = OUTPUT_DAT + "/01_emfs/emf_list_str"
    output:
        OUTPUT_REP + "/emf_summary.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "emf_summary.Rmd"

#-------------------------------------------------------------------------------
"""
# Re-clustering data with only the genes that are 
#
# - non-differentially expressed 
# - conserved marker genes
# - core genes (rename later)
#
#To control that the genes are useful
"""

rule clustering_ctrl:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        core_cons_list = rules.export_emf_genes.output
    params: 
        k_louvain = config["values"]["02_sce_anno"]["k_louvain"]
    output:
        sce_output = OUTPUT_DAT + "/02_clst/sce_{fraction}"
    script:
        "scripts/02_clustering_ctrl.R"    
        
rule clustering_refs:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        sce_baccin_in = config["base"] + config["metadata_paths"]["ref_baccin_sce"],
        sce_dolgalev = config["base"] + config["metadata_paths"]["ref_dolgalev_sce"],
        sce_dahlin = config["base"] + config["metadata_paths"]["ref_dahlin_sce"],
        core_cons_list = rules.export_emf_genes.output
    params: 
        k_louvain = config["values"]["02_sce_anno"]["k_louvain_refs"]
    output:
        sce_baccin_out = OUTPUT_DAT + "/03_cref/sce_baccin_{fraction}",
        sce_ref2_out = OUTPUT_DAT + "/03_cref/sce_ref2_{fraction}"
    script:
        "scripts/03_clustering_refs.R" 

rule clustering_ctrl_report:
    input: 
        sce_input = rules.clustering_ctrl.output,
        sce_baccin = rules.clustering_refs.output.sce_baccin_out,
        sce_ref2 = rules.clustering_refs.output.sce_ref2_out
    output:
        OUTPUT_REP + "/clustering_ctrl_report_{fraction}.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "clustering_ctrl_report.Rmd"

#-------------------------------------------------------------------------------

