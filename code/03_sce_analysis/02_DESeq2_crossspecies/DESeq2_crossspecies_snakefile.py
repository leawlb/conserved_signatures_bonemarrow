#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/02_DESeq2_crossspecies"

COLORS = config["base"] + config["metadata_paths"]["colors"]

CELL_TYPES_EXCLUDE = config["values"]["03_sce_analysis"]["cell_types_exclude"]
print(CELL_TYPES_EXCLUDE)
SV_PATH = config["base"] + config["metadata_paths"]["sources_variation"]["analysis_species"]
print(SV_PATH)

METADATA = pd.read_csv(config["table"])
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

# make cell types list from annotation
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
print(celltypes)
celltypes_hsc = celltypes[0:12]
celltypes_str = celltypes[12:22]

print(celltypes_hsc)
print(celltypes_str)

PADJ_CUTOFF = config["values"]["02_sce_anno"]["ndge_padj_cutoff"]
FC_CUTOFF = config["values"]["02_sce_anno"]["ndge_fc_cutoff"]
tf = "PC_" + str(PADJ_CUTOFF) + "_FC_" + str(FC_CUTOFF)
print(tf)

#-------------------------------------------------------------------------------

targets = []

for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_desq/deseq_" + f]
  targets = targets + [OUTPUT_DAT + "/02_dsqc/rlog_" + f]
  targets = targets + [OUTPUT_DAT + "/02_dsqc/sva_" + f]
  targets = targets + [OUTPUT_DAT + "/04_tdsq/deseq_" + f]
  targets = targets + [OUTPUT_DAT + "/05_nres/" + tf + "/res_" + f + "_celltype"]
  targets = targets + [OUTPUT_DAT + "/05_nres/" + tf + "/res_" + f + "_celltype_dfs"]
  targets = targets + [OUTPUT_DAT + "/05_nres/" + tf + "/shared_genes_" + f + "_celltypes"]

for c in celltypes_hsc:
  targets = targets + [OUTPUT_DAT + "/03_sepd/hsc_" + c + "-sep"]
  targets = targets + [OUTPUT_REP + "/sva/sva_report_hsc_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_hsc_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_hsc_" + c + ".html"] 
# 
for c in celltypes_str:
  targets = targets + [OUTPUT_DAT + "/03_sepd/str_" + c + "-sep"]
  targets = targets + [OUTPUT_REP + "/sva/sva_report_str_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_str_" + c + ".html"] 
  targets = targets + [OUTPUT_REP + "/ndge/ndge_report_str_" + c + ".html"] 

#-------------------------------------------------------------------------------

wildcard_constraints: 
  fraction="[a-z]+"

localrules: all  

rule all: 
  input:
      targets

#-------------------------------------------------------------------------------
"""
# DESeq pre-processing and QC + Reports
"""

"""
# aggregate pseudobulks per cell type per sample and convert to DESeq2 object
# output = list of DESeq2 objects generated per fraction
# one DESeq2 object/list item for each cell type
"""
rule aggregate_convert:
  input:
      sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10"
  resources:
      mem_mb=30000,
      queue="medium-debian"
  threads: 4
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
    resources:
        mem_mb=5000,
        queue="medium-debian"
    threads: 4
    output:
        rlog = OUTPUT_DAT + "/02_dsqc/rlog_{fraction}",
        sva = OUTPUT_DAT + "/02_dsqc/sva_{fraction}"
    script:
        "scripts/02_prep_qc_deseq.R"  
        
# dummy rule to allow easy use of cell types as wildcard for reports
output = expand(OUTPUT_DAT + "/03_sepd/hsc_{celltype}-sep", celltype = celltypes_hsc)
output = output + expand(OUTPUT_DAT + "/03_sepd/str_{celltype}-sep", celltype = celltypes_str)
rule separate_sce:
    input: 
        sce_input = expand(OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10", fraction = fractions)
    resources:
        mem_mb=50000,
        queue="medium-debian"
    threads: 1
    output:
        output = output
    script:
        "scripts/03_separate_dummy.R"
        
#-------------------------------------------------------------------------------

# check hidden sources of variations and decide which ones to add to design
rule celltype_sva_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_{celltype}-sep",
        dsq_list = rules.aggregate_convert.output,
        sva_list = rules.qc_deseq.output.sva
    resources:
        mem_mb=5000,
        queue="medium-debian"
    threads: 4
    output:
        OUTPUT_REP + "/sva/sva_report_{fraction}_{celltype}.html"
    params:
        plotting = "../../source/plotting.R"
    script:
        "sva_report.Rmd"
        
# check quality of the pseudo-bulks at cell type level (annotated)       
rule celltype_bulk_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_{celltype}-sep",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        rlog = rules.qc_deseq.output.rlog
    resources:
        mem_mb=25000,
        queue="medium-debian"
    threads: 4
    output:
        OUTPUT_REP + "/bulk/bulk_quality_report_{fraction}_{celltype}.html"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R",
        cts_exclude = CELL_TYPES_EXCLUDE
    script:
        "bulk_quality_report.Rmd"

#-------------------------------------------------------------------------------



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
        sva = rules.qc_deseq.output.sva,
        sv_path = SV_PATH
    resources:
        mem_mb=5000,
        queue="medium-debian"
    threads: 4
    output:
        deseq_output = OUTPUT_DAT + "/04_tdsq/deseq_{fraction}"
    script:
        "scripts/04_deseq_species.R"    
        
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
        species = species,
        cts_exclude = CELL_TYPES_EXCLUDE
    resources:
        mem_mb=5000,
        queue="medium-debian"
    threads: 4
    output:
        celltype_res = OUTPUT_DAT + "/05_nres/" + tf + "/res_{fraction}_celltype",
        celltype_resdf_list = OUTPUT_DAT + "/05_nres/" + tf + "/res_{fraction}_celltype_dfs",
        celltype_shared_genes_list = OUTPUT_DAT + "/05_nres/" + tf + "/shared_genes_{fraction}_celltypes"
    script:
        "scripts/05_export_ndge_results.R"    


#-------------------------------------------------------------------------------

"""   
# Make Reports
# This part of the script was adjusted from 02_sce_anno/02_nDGE/
"""

rule celltype_ndge_report:
    input: 
        sep = OUTPUT_DAT + "/03_sepd/{fraction}_{celltype}-sep",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        celltype_res = rules.export_results_ndge.output.celltype_res,
        celltype_shared_genes_list = rules.export_results_ndge.output.celltype_shared_genes_list,
    output:
        OUTPUT_REP + "/ndge/ndge_report_{fraction}_{celltype}.html"
    resources:
        mem_mb=25000,
        queue="medium-debian"
    threads: 4
    params:
        padj_cutoff = PADJ_CUTOFF,
        fc_cutoff = FC_CUTOFF,
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R",
        cts_exclude = CELL_TYPES_EXCLUDE
    script:
        "ndge_celltype_report.Rmd"

