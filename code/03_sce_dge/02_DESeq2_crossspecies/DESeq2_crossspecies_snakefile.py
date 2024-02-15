#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

# paths from config
OUTPUT_BASE = config["base"] + config["data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies"
OUTPUT_REP = OUTPUT_BASE + "/reports/03_sce_analysis/02_DESeq2_crossspecies"

COLORS = config["base"] + config["metadata_paths"]["colors"]

# objects from config
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
print(fractions)

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
#celltypes_hsc = celltypes[0:1]
#celltypes_str = celltypes[12:13]

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
  targets = targets + [OUTPUT_DAT + "/06_sepd/hsc_" + c + "-sep"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_hsc_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/dge/dge_report_hsc_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_hsc_" + c + ".html"] 

for c in celltypes_str:
  targets = targets + [OUTPUT_DAT + "/06_sepd/str_" + c + "-sep"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_str_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/dge/dge_report_str_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_str_" + c + ".html"] 


for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_desq/deseq_" + f]
  targets = targets + [OUTPUT_DAT + "/02_dsqc/rlog_" + f]
  targets = targets + [OUTPUT_DAT + "/02_dsqc/sva_" + f]
  targets = targets + [OUTPUT_DAT + "/03_tdsq/deseq_" + f]
  targets = targets + [OUTPUT_DAT + "/04_dres/res_" + f + "_celltype"]
  targets = targets + [OUTPUT_DAT + "/04_dres/res_" + f + "_celltype_dfs"]
  targets = targets + [OUTPUT_DAT + "/05_nres/" + tf + "/res_" + f + "_species"]
  targets = targets + [OUTPUT_DAT + "/05_nres/" + tf + "/res_" + f + "_celltype"]
  targets = targets + [OUTPUT_DAT + "/05_nres/" + tf + "/res_" + f + "_celltype_dfs"]
  targets = targets + [OUTPUT_DAT + "/05_nres/" + tf + "/res_" + f + "_celltype_shared"]

#print(targets)

#-------------------------------------------------------------------------------

wildcard_constraints: #should not contain ".", see rule cellranger_count
  fraction="[a-z]+"

localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules
  input:
      targets

#-------------------------------------------------------------------------------
"""
# DESeq pre-processing and QC 
"""

"""
# aggregate pseudobulks per cell type per sample and convert to DESeq2 object
# output = list of DESeq2 objects generated per fraction
# one DESeq2 object/list item for each cell type
"""
rule aggregate_convert:
  input:
      sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10"
  output:
      deseq_output = OUTPUT_DAT + "/01_desq/deseq_{fraction}"
  script:
      "scripts/01_aggregate_convert.R"


"""
# QC of aggregated DESeq2 objects
# produce rlog-transformed values purely for QC visualisation
# identify hidden sources of variations for batch correction of DESeq2 objects
"""
rule qc_deseq:
    input:
        deseq_input = rules.aggregate_convert.output
    output:
        rlog = OUTPUT_DAT + "/02_dsqc/rlog_{fraction}",
        sva = OUTPUT_DAT + "/02_dsqc/sva_{fraction}"
    script:
        "scripts/02_prep_qc_deseq.R"  

#-------------------------------------------------------------------------------

"""
# Cross-species analysis
"""

"""
# preprocessing/transformation/normalisation of DESeq2 Objects for 
# cross-species DGE and nDGE analysis
# SVs are taken into account for desiggn, as batch correction
"""
rule preprocessing_deseq:
    input:
        deseq_input = rules.aggregate_convert.output,
        sva = rules.qc_deseq.output.sva
    output:
        deseq_output = OUTPUT_DAT + "/03_tdsq/deseq_{fraction}"
    script:
        "scripts/03_deseq_species.R"    

"""
# exporting DGE results
# per cell type, comparing across species
# using classic DESeq2 function for each pairwise comparison across species
# for each cell type
"""
rule export_results_dge:
    input:
        deseq_input = rules.preprocessing_deseq.output
    params:
        species = species
    output:
        celltype_res = OUTPUT_DAT + "/04_dres/res_{fraction}_celltype",
        celltype_res_dfs = OUTPUT_DAT + "/04_dres/res_{fraction}_celltype_dfs",
    script:
        "scripts/04_export_dge_results_species.R"    
        
"""
# exporting nDGE results 
# using alternative hypothesis with stated logFC and pval cut-offs
# exporting lists of ndges that are shared across specied
# ndges will be used in downstream analysis
"""
rule export_results_ndge:
    input:
        deseq_input = rules.preprocessing_deseq.output
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        species = species
    output:
        species_res = OUTPUT_DAT + "/05_nres/" + tf + "/res_{fraction}_species",
        celltype_res = OUTPUT_DAT + "/05_nres/" + tf + "/res_{fraction}_celltype",
        celltype_res_dfs = OUTPUT_DAT + "/05_nres/" + tf + "/res_{fraction}_celltype_dfs",
        celltype_res_list_shared = OUTPUT_DAT + "/05_nres/" + tf + "/res_{fraction}_celltype_shared"
    script:
        "scripts/05_export_ndge_results.R"    

#-------------------------------------------------------------------------------
"""
Dummy rule to allow easy use of cell types as wildcard for reports
Unfortunately this takes quite long
"""

output = expand(OUTPUT_DAT + "/06_sepd/hsc_{celltype}-sep", celltype = celltypes_hsc)
output = output + expand(OUTPUT_DAT + "/06_sepd/str_{celltype}-sep", celltype = celltypes_str)
rule separate_sce:
    input: 
        sce_input = expand(OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10", fraction = fractions)
    output:
        output = output
    script:
        "scripts/06_separate_dummy.R"
        
#-------------------------------------------------------------------------------

"""   
# Make Reports
# This part of the script was adjusted from 02_sce_anno/02_nDGE/
"""

rule celltype_bulk_report:
    input: 
        sep = OUTPUT_DAT + "/06_sepd/{fraction}_{celltype}-sep",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        rlog = rules.qc_deseq.output.rlog
    output:
        OUTPUT_REP + "/bulk/bulk_quality_report_{fraction}_{celltype}.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "bulk_quality_report.Rmd"
 

rule celltype_dge_report:
    input: 
        sep = OUTPUT_DAT + "/06_sepd/{fraction}_{celltype}-sep",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        celltype_res = rules.export_results_dge.output.celltype_res,
    output:
        OUTPUT_REP + "/dge/dge_report_{fraction}_{celltype}.html"
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "dge_celltype_report.Rmd"


rule celltype_ndge_report:
    input: 
        sep = OUTPUT_DAT + "/06_sepd/{fraction}_{celltype}-sep",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        celltype_res_list_shared = rules.export_results_ndge.output.celltype_res_list_shared,
        celltype_res_list_shared_orig = OUTPUT_BASE + "/sce_objects/02_sce_anno/08_dres/PC_0.05_FC_1.5/res_{fraction}_cluster_shared",
        celltype_res = rules.export_results_ndge.output.celltype_res,
        celltype_res_orig = OUTPUT_BASE + "/sce_objects/02_sce_anno/08_dres/PC_0.05_FC_1.5/res_{fraction}_cluster",
    output:
        OUTPUT_REP + "/ndge/ndge_report_{fraction}_{celltype}.html"
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "ndge_celltype_report.Rmd"
