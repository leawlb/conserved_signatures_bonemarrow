#!/bin/python 

# Monster Snakefile of entire CCI process to compare conditions.
# For more details, check out 04_cci_prep and 06_cci_analysis.

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/02_cci_alt/02_reverse"
OUTPUT_REP = OUTPUT_BASE + "/cci_objects/reports/02_cci_alt/02_reverse"

COLORS = config["base"] + config["metadata_paths"]["colors"]

LRDB_OUT =  config["base"] + config["metadata_paths"]["lrdb_out"]
ASSIGNMENT = config["base"] + config["metadata_paths"]["assignment_reverse"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
age = get_list(metadata = METADATA, column = "Age_ID")

VALUES = config["values"]["05_cci_prep"]


#-------------------------------------------------------------------------------

targets = []
for s in species:
  for a in age:
    targets = targets + [OUTPUT_DAT + "/01_asnm_rev/sce_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/02_subs_rev/sce_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/03_down_rev/sce_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/04_calc_rev/cci_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/05_mtrc_rev/cci_metrics_" + s + "_" + a]
    targets = targets + [OUTPUT_DAT + "/06_mrge_rev/cci_all"]
    targets = targets + [OUTPUT_DAT + "/07_ed_rev/cci_dist_all"]

targets = targets + [OUTPUT_REP + "/reverse_summary.html"]
    
#-------------------------------------------------------------------------------

localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

"""
Prepare SCE objects for CCI using the reversed assignment.
Subsetting and Downsampling steps should actually be identical
"""
        
rule celltype_assignment_rev:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/01_sprp/sce_{species}_{age}-01",
        assignment = ASSIGNMENT
    params:
        min_cells = VALUES["min_cells"]
    output:
        sce_output = OUTPUT_DAT + "/01_asnm_rev/sce_{species}_{age}"
    script:
        "scripts/01_add_assignment_rev.R" 
        
rule subset_sce_rev:
    input:
        sce_input = rules.celltype_assignment_rev.output,
        lrdb_input = LRDB_OUT
    output:
        sce_output = OUTPUT_DAT + "/02_subs_rev/sce_{species}_{age}",
    script:
        "scripts/02_subset_sce_rev.R" 
        
rule downsample:
    input:
        sce_input_path = expand(rules.subset_sce_rev.output, species = species, age = age),
        assignment = ASSIGNMENT
    output:
        sce_output_path = expand(OUTPUT_DAT + "/03_down_rev/sce_{species}_{age}", species = species, age = age)
    script:
        "scripts/03_downsample_ct_rev.R" 
        
"""

Calculate CCIs under main conditions, but with reversed assignments.
To test if generally observed trends on the impact of ageing/species are
true also for this direction of signalling

The following rules summarize the entire process of CCI calculation.

"""

rule calculate_rev:
    input:
        sce_input = OUTPUT_DAT + "/03_down_rev/sce_{species}_{age}",
        lrdb = LRDB_OUT
    output:
        cci_output = OUTPUT_DAT + "/04_calc_rev/cci_{species}_{age}"
    params:
        main_functions = "../../source/cci_functions_calculation.R",
        alternative_functions = "../../source/cci_functions_alt_conds.R",
        min_perc = VALUES["min_perc"],
        top_level = VALUES["top_level"],
        assay_use = "downsampled",
        condition = "all preprocessing steps"
    script:
        "scripts/04_calculate_rev.R" 

"""

Calculate CCI metrics for reversed assignments
The following rules summarizes the entire process of CCI metrics calculation.

"""

rule metrics_rev:
    input:
        sce_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_{age}-05",
        cci_input =  OUTPUT_DAT + "/04_calc_rev/cci_{species}_{age}"
    output:
        metrics_lists = OUTPUT_DAT + "/05_mtrc_rev/cci_metrics_{species}_{age}"
    params:
        metrics_functions = "../../source/cci_functions_metrics.R"
    script:
        "scripts/05_metrics_rev.R" 
        
"""
Calculate Euclidian distances for reversed assignments after merging
"""

rule merge_cci_rev:
    input: 
        cci_input_path = expand(OUTPUT_DAT + "/04_calc_rev/cci_{s}_{a}", a = age, s = species)
    output:
        cci_output = OUTPUT_DAT + "/06_mrge_rev/cci_all"
    script:
        "scripts/06_merge_rev.R"   

rule calculate_ed_rev:
    input: 
        cci_input = OUTPUT_DAT + "/06_mrge_rev/cci_all"
    output:
        cci_dist_output = OUTPUT_DAT + "/07_ed_rev/cci_dist_all"
    script:
        "scripts/07_ed_rev.R"   

"""
Reports
"""

rule rev_report:
    input:
        metrics_path = expand(rules.metrics_rev.output.metrics_lists, age = age, species = species),
        cci_dist_all = rules.calculate_ed_rev.output
    output:
        OUTPUT_REP + "/reverse_summary.html"
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"    
    script:
        "cci_rev_assignment_report.Rmd" 
