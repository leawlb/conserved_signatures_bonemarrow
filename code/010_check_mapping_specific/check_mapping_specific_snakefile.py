#!/bin/python 

#-------------------------------------------------------------------------------

import pandas as pd

# paths from config
METADATA_PATH = config["metadata"]["table"]
OUTPUT_BASE_PATH = config["paths"]["output_dir"]
TABLES_PATH = config["metadata"]["color_tables"]

# objects from config
IDENTIFIERS = config["metadata"]["identifiers"]
METADATA = pd.read_csv(config["metadata"]["table"])
VALUES =  config["values"]["preprocessing"]

def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

species = get_list(metadata = METADATA, column = "Species_ID")
fractions = get_list(metadata = METADATA, column = "Fraction_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

fractions = ["hsc"]

# collect all paths for all possible outputs/targets, required for rule all
targets = []
for s in species:
  for f in fractions:
    targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/14_check_mapping_specific/sce_tog_" + s + "_" + f + "-14"]
    targets = targets + [OUTPUT_BASE_PATH + "/sce_objects/14_check_mapping_specific/og_markergenes_" + s + "_" + f]
    targets = targets + [OUTPUT_BASE_PATH + "/reports/010_check_mapping_specific/mapping_species_report_" + s + "_" + f + ".html"]
print(targets)
# local execution of non-demanding rules
localrules: all  

#-------------------------------------------------------------------------------

# define rules
rule all: # must contain all possible output paths from all rules below
    input:
        targets

#-------------------------------------------------------------------------------
# prepare fg (four genomes) and og (one genome) objects for comparison
rule prepare_objects:
    input:
        sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_species/sce_{species}_{fraction}-07",
        sce_fg_path = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/fourgenomes/sce_objects/01_cellranger_output/"
    output:
        sce_tog = OUTPUT_BASE_PATH + "/sce_objects/14_check_mapping_specific/sce_tog_{species}_{fraction}-14"
    params:
        nr_hvgs = VALUES["nr_hvgs"],
        individuals = individuals,
    script:
        "scripts/14_prepare_objects.R"
    
rule get_mapping_diff_markergenes:
    input:
        sce_tog = rules.prepare_objects.output
    output:
        markergenes_mapping = OUTPUT_BASE_PATH + "/sce_objects/14_check_mapping_specific/og_markergenes_{species}_{fraction}"
    params:
        nr_hvgs = VALUES["nr_hvgs"],
        individuals = individuals,
    script:
        "scripts/14_markergenes_mapping.R"

# construct report files to monitor QC on sample level

rule make_species_reports:
    input:
        sce_tog = rules.prepare_objects.output,
        markergenes_mapping = rules.get_mapping_diff_markergenes.output,
        # TODO make stromal gene list
        gene_list = config["metadata"]["gene_list_hsc"],
        dge_results = OUTPUT_BASE_PATH + "/sce_objects/13_DGE_analysis_flayer/PC_0.001_FC_1.3/DGE/res_{fraction}_{species}_specific_df-11.csv",
    output:
        OUTPUT_BASE_PATH + "/reports/010_check_mapping_specific/mapping_species_report_{species}_{fraction}.html"
    params:
        nr_hvgs = VALUES["nr_hvgs"]
    script:
        # .Rmd files should be stored in snakefile working directory
        "check_mapping_specific.Rmd" 
