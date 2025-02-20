#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------
# CONFIG
# get all paths and objects from config

# base path
OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
print(OUTPUT_BASE)

# directory for storing data
OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/01_sce_prep"
# directory for storing reports
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/01_sce_prep/01_preprocessing"

# directories of raw input data
CELLRANGER_OUT = config["cellranger_output"] 
FOURGENOMES_OUT = config["fourgenomes_output"]

# load metadata table to obtain species, individuals, fractions to be used
METADATA = pd.read_csv(config["table"])

# function to extract the required variables
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

# get required variables for targets and wildcards
species = get_list(metadata = METADATA, column = "Species_ID")
individuals = get_list(metadata = METADATA, column = "Object_ID")

# parameter values stored in config
VALUES = config["values"]["01_sce_prep"]

#-------------------------------------------------------------------------------

# collect all paths for outputs/targets, required for rule all
targets = []
for s in species:
  for i in individuals:
    if s in i:
      targets = targets + [OUTPUT_DAT + "/01_drop/" + s + "/sce_" + i + "-01"]
      targets = targets + [OUTPUT_DAT + "/02_mapp/" + s + "/dmgs_" + i]
      targets = targets + [OUTPUT_DAT + "/04_outl/" + s + "/sce_" + i + "-04"]
      targets = targets + [OUTPUT_DAT + "/05_norm/" + s + "/sce_" + i + "-05"]
      targets = targets + [OUTPUT_DAT + "/06_dimr/" + s + "/sce_" + i + "-06"]
      targets = targets + [OUTPUT_REP + "/qc/" + s + "/preprocessing_qc_report_" + i + ".html"]
      targets = targets + [OUTPUT_REP + "/dmgs/" + s + "/preprocessing_dmg_report_" + i + ".html"]

targets = targets + [OUTPUT_DAT + "/03_dmgs/dmgs_list_all"]
      
if config["run_preprocessing_summary"]:
  targets = targets + [OUTPUT_REP + "/qc/preprocessing_qc_summary.html"]

#-------------------------------------------------------------------------------

# local execution of non-demanding rules
localrules: all  

# define rules
rule all: # must contain all possible output paths from all rules below and must always be the first rule
    input:
        targets
        
#-------------------------------------------------------------------------------
# remove empty droplets, doublets, duplicated genes from main data 
rule remove_droplets:
    input: 
        sce_input = CELLRANGER_OUT + "/{species}/sce_{individual}-01"
    output:
        sce_output = OUTPUT_DAT + "/01_drop/{species}/sce_{individual}-01"
    resources:
        mem_mb=10000,
        queue="medium"
    threads: 2
    params:
        cutoff_umis = VALUES["cutoff_umis"],
        cutoff_doublets = VALUES["cutoff_doublets"]
    script:
        "scripts/01_remove_droplets.R"

"""
# Get the differentially mapped genes (dmgs) for each sample by comparing
# the expression levels between the main data, and data annotated with species-
# specific genomes

# ensembl_list_mspr contains gene annotations for mspr-specific genome-
# annotated SCE objects
"""
rule get_sample_dmgs:
    input:
        sce_input = rules.remove_droplets.output,
        sce_fg = FOURGENOMES_OUT +  "/{species}/sce_{individual}-01",
        ensembl_list_mspr = config["ensembl_mspr"]
    output:
        dmgs = OUTPUT_DAT + "/02_mapp/{species}/dmgs_{individual}"
    resources:
        mem_mb=15000,
        queue="short"
    threads: 2
    params:
        nr_hvgs = config["values"]["nr_hvgs"],
        logFC_sample_dmgs = VALUES["logFC_sample_dmgs"],
        minPC_sample_dmgs = VALUES["minPC_sample_dmgs"]
    script:
        "scripts/02_sample_dmgs.R" 
    
# merge all dmgs into one list for exclusion 
merge_dmgs_input = []
for s in species:
  for i in individuals:
    if s in i:
      merge_dmgs_input = merge_dmgs_input + expand(rules.get_sample_dmgs.output, species = s, individual = i)

rule merge_dmgs:
    input:
        dmgs = merge_dmgs_input
    resources:
        mem_mb=100,
        queue="short"
    threads: 1
    output:
        dmg_list = OUTPUT_DAT + "/03_dmgs/dmgs_list_all"
    script:
        "scripts/03_merge_dmgs.R"
        
# remove cells with low quality, also remove dmgs 
rule remove_outliers_dmgs:
    input: 
        sce_input = rules.remove_droplets.output,
        dmg_list = rules.merge_dmgs.output
    output:
        sce_output = OUTPUT_DAT + "/04_outl/{species}/sce_{individual}-04"
    resources:
        mem_mb=2000,
        queue="medium"
    threads: 2
    params:
        cutoff_sum = VALUES["cutoff_sum"],
        cutoff_detected = VALUES["cutoff_detected"],
        cutoff_mitos = VALUES["cutoff_mitos"]
    script:
        "scripts/04_remove_outliers_dmgs.R"
        
# normalize expression on sample level
rule normalize_expr:
    input: 
        sce_input = rules.remove_outliers_dmgs.output
    output:
        sce_output = OUTPUT_DAT + "/05_norm/{species}/sce_{individual}-05"
    resources:
        mem_mb=2000,
        queue="short"
    threads: 2
    script:
        "scripts/05_normalize_expr.R"
        
# extract hvgs and reduce dimensions for QC visualisation on sample level 
# use own function that summarizes dimensionaility reduction steps
rule reduce_dims:
    input: 
        sce_input = rules.normalize_expr.output
    output:
        sce_output = OUTPUT_DAT + "/06_dimr/{species}/sce_{individual}-06" 
    resources:
        mem_mb=2000,
        queue="short"
    threads: 2
    params:
        nr_hvgs = config["values"]["nr_hvgs"],
        functions = "../../source/sce_functions.R"
    script:
        "scripts/06_reduce_dims.R"
        
#-------------------------------------------------------------------------------
# REPORTS

# report file for dmgs 
rule make_dmg_reports:
    input:
        sce_input = rules.remove_droplets.output,
        sce_fg = FOURGENOMES_OUT + "/{species}/sce_{individual}-01",
        dmgs = rules.get_sample_dmgs.output,
        dmg_list = rules.merge_dmgs.output,
        ensembl_list_mspr = config["ensembl_mspr"]
    output:
        OUTPUT_REP + "/dmgs/{species}/preprocessing_dmg_report_{individual}.html"
    resources:
        mem_mb=20000,
        queue="medium"
    threads: 2
    params:
        nr_hvgs = config["values"]["nr_hvgs"],
        plotting = "../../source/plotting.R" # path to source file
    script:
        # .Rmd files should be stored in snakefile working directory, 
        # which is the directory of the snakefile
        # the working directory of scripts in the /scripts folder 
        # is also the directory of the snakefile
        "preprocessing_dmg_reports.Rmd" 

# report file for QC
rule make_qc_reports:
    input:
        sce_input = CELLRANGER_OUT + "/{species}/sce_{individual}-01",
        sce_drop = rules.remove_droplets.output,
        sce_outl = rules.remove_outliers_dmgs.output,
        sce_norm = rules.normalize_expr.output,
        sce_dimr = rules.reduce_dims.output
    output:
        OUTPUT_REP + "/qc/{species}/preprocessing_qc_report_{individual}.html"
    resources:
        mem_mb=20000,
        queue="medium"
    threads: 2
    params:
        cutoff_umis = VALUES["cutoff_umis"],
        cutoff_doublets = VALUES["cutoff_doublets"],
        cutoff_sum = VALUES["cutoff_sum"],
        cutoff_detected = VALUES["cutoff_detected"],
        cutoff_mitos = VALUES["cutoff_mitos"],
        nr_hvgs = config["values"]["nr_hvgs"],
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R" 
    script:
        "preprocessing_qc_reports.Rmd" 
        
"""
Make one summary report on all files
Because this rule is dependent on previous outputs but this is not reflected
in input, it's generally turned off in config to avoid premature job submission
"""
print(config["run_preprocessing_summary"])
if config["run_preprocessing_summary"]:
  rule make_summary_report:
      input:
          sce_input_path = OUTPUT_DAT + "/01_drop/"
      output:
          OUTPUT_REP + "/qc/preprocessing_qc_summary.html"
      resources:
          mem_mb=20000,
          queue="long"
      threads: 2
      params:
          cutoff_sum = VALUES["cutoff_sum"],
          cutoff_detected = VALUES["cutoff_detected"],
          cutoff_mitos = VALUES["cutoff_mitos"],
          individuals = individuals,
          metadata = config["table"],
          plotting = "../../source/plotting.R",
          species = species
      script:
          "preprocessing_qc_summary.Rmd" 
