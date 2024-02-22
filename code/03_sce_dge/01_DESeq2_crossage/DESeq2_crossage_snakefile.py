#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

# paths from config
OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/01_DESeq2_crossage"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/01_DESeq2_crossage"

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

#-------------------------------------------------------------------------------

# construct paths for all possible outputs/targets, required for rule all
targets = []

for f in fractions:
  for s in species:
    targets = targets + [OUTPUT_DAT + "/01_desq/deseq_" + s + "_" + f]
    targets = targets + [OUTPUT_DAT + "/02_dsqc/rlog_" + s + "_" + f]
    targets = targets + [OUTPUT_DAT + "/02_dsqc/sva_" + s + "_"+ f]
    targets = targets + [OUTPUT_DAT + "/03_tdsq/deseq_" + s + "_" + f]
    targets = targets + [OUTPUT_DAT + "/04_dres/res_" + s + "_" + f]
    targets = targets + [OUTPUT_REP + "/bulk/bulk_quality_report_" + s + "_" + f + ".html"]
    targets = targets + [OUTPUT_REP + "/dge/dge_report_" + s + "_" + f + ".html"]

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
# Cross-age analysis
"""

#-------------------------------------------------------------------------------

"""
# DESeq pre-processing and QC 
"""

"""
# aggregate pseudobulks per cell type and convert to DESeq2 object
# output = list of DESeq2 objects generated per fraction
# one DESeq2 object/list item for each cell type
# one list for each species
"""
rule aggregate_convert:
  input:
      sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10"
  output:
      deseq_output = expand(OUTPUT_DAT + "/01_desq/deseq_{species}_{{fraction}}", species = species)
  script:
      "scripts/01_aggregate_convert.R"

"""
# QC of aggregated DESeq2 objects
# produce rlog-transformed values purely for QC visualisation
# identify hidden sources of variations for batch correction of DESeq2 objects
"""
rule qc_deseq:
    input:
        deseq_input = OUTPUT_DAT + "/01_desq/deseq_{species}_{fraction}"
    output:
        rlog = OUTPUT_DAT + "/02_dsqc/rlog_{species}_{fraction}",
        sva = OUTPUT_DAT + "/02_dsqc/sva_{species}_{fraction}"
    script:
        "scripts/02_prep_qc_deseq.R"  

#-------------------------------------------------------------------------------

"""
# preprocessing/transformation/normalisation of DESeq2 Objects for 
# cross-species DGE and nDGE analysis
# SVs are taken into account for desiggn, as batch correction
"""
rule preprocessing_deseq:
    input:
        deseq_input = OUTPUT_DAT + "/01_desq/deseq_{species}_{fraction}",
        sva = rules.qc_deseq.output.sva
    output:
        deseq_output = OUTPUT_DAT + "/03_tdsq/deseq_{species}_{fraction}"
    script:
        "scripts/03_deseq_ages.R"    

"""
# exporting DGE results
# per cell type, comparing across species
# using classic DESeq2 function for each pairwise comparison across species
# for each cell type
"""
rule export_results_dge:
    input:
        deseq_input = rules.preprocessing_deseq.output
    output:
        res_list = OUTPUT_DAT + "/04_dres/res_{species}_{fraction}",
    script:
        "scripts/04_export_dge_results.R"    

#-------------------------------------------------------------------------------

"""   
# Make Reports
"""

rule species_bulk_report:
    input: 
        rlog = rules.qc_deseq.output.rlog
    output:
        OUTPUT_REP + "/bulk/bulk_quality_report_{species}_{fraction}.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "species_bulk_quality_report.Rmd"

rule dge_report:
    input: 
        res_list = rules.export_results_dge.output.res_list,
    output:
        OUTPUT_REP + "/dge/dge_report_{species}_{fraction}.html"
    params:
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "dge_report.Rmd"
