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

# collect all paths for all possible outputs/targets, required for rule all
targets = []
for s in species:
  for ind in individuals:
      if s in ind:
          targets = targets + [OUTPUT_BASE_PATH + "/reports/010_check_mapping_specific/" + s + "/" + f + "/mapping_species_report.html" for f in fractions]

print(targets)
# local execution of non-demanding rules
localrules: all  

#-------------------------------------------------------------------------------

# define rules
rule all: # must contain all possible output paths from all rules below
    input:
        targets

#-------------------------------------------------------------------------------
# construct report files to monitor QC on sample level

rule make_species_reports:
    input:
        sce_07 = OUTPUT_BASE_PATH + "/sce_objects/07_mrge_species/sce_{species}_{fraction}-07",
        sce_fg_path = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/fourgenomes/sce_objects/01_cellranger_output/"
    output:
        OUTPUT_BASE_PATH + "/reports/010_check_mapping_specific/{species}/{fraction}/mapping_species_report.html"
    params:
        nr_hvgs = VALUES["nr_hvgs"],
        color_tables = TABLES_PATH,
        individuals = individuals,
        gene_list_hsc_all = "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/gene_list_hsc_all.csv"
    script:
        # .Rmd files should be stored in snakefile working directory
        "check_mapping_specific.Rmd" 
