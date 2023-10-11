#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["data_paths"]["main"]
OUTPUT_DAT = OUTPUT_BASE + "/cci_objects/03_cci_analysis"
OUTPUT_REP = OUTPUT_BASE + "/reports/06_cci_analysis/03_cci_permutation"

COLORS = config["base"] + config["metadata_paths"]["colors"]

METADATA = pd.read_csv(config["base"] + config["metadata_paths"]["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)
  
species = get_list(metadata = METADATA, column = "Species_ID")
ages = get_list(metadata = METADATA, column = "Age_ID")

VALUES = config["values"]["04_cci_prep"]
LRDB = config["base"] + config["metadata_paths"]["lrdb_out"]

#-------------------------------------------------------------------------------

targets = []

for s in species:
  targets = targets + [OUTPUT_DAT + "/04_pres_age/res_lists_ages_" + s]
  targets = targets + [OUTPUT_REP + "/report_perm_age/" + s + "_report.html"]
  
  for a in ages:
    targets = targets + [OUTPUT_DAT + "/03_perm_age/perm_scores_" + s + "_" + a]

targets = targets + [OUTPUT_DAT + "/05_psum_age/summary_df"]
targets = targets + [OUTPUT_REP + "/report_perm_age/perm_summary.html"]


#-------------------------------------------------------------------------------


localrules: all  

rule all: 
    input:
        targets
 
#-------------------------------------------------------------------------------

"""
Permute age labels within species

This scripts summarizes the entire process of generating CCI objects from down-
sampled SCE objects as shown in 04_cci_prep/02_cci_calculation.
The main functions for CCI calculation were adjusted for scRNAseq data
from Adrien Jolly's CellInteractionScores Method:
https://github.com/AdrienJolly/CellInteractionScores

1. Permutation
Age labels in SCE objects are permuted and CCI objects are calculated 
x n = iterations for each permuted age.
All n calculated scores per interacting  cell type pair (ctp) per age are
extracted and returned in a list of data frames.

2. Testing
Score deltas (score old - score yng) are calculated to obtain the test statistic 
and permutation distribution per interaction per cell type pair.
P values are estimated for each interaction in each cell type pair.
Permutation distribution and other metrics are visualised in reports. 

It is important to note that pvals are only estimated.
Using sets of permutations instead of all possible permutations results in
estimated p values that are generally understated, and also results in more
tests that reject the Null Hypothesis. 
It also results in p values that are exactly 0 which cannot be adjusted for 
multiple testing.
But using estimated p values helps identify interactions, cell types, or
cell type pairs of interest.

Phipson, B., and Smyth, G. K., Stat. Appl. Genet. Molec. Biol. (2010)

"""

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""

Permute age labels and calculate CCI objects

Age labels from young and old SCE objects are permuted (n times).
Then for each permuted SCE object ("yng" or "old"), a CCI object is calculated
as described in 04_cci_prep/02_cci_calculation.
Next, permuted "yng" and "old" scores are extracted and stored in dataframes 
with rows = interactions and cols = iterations.
One dataframe for each interacting cell type pair.

Output = list of dataframes containing permuted scores, one list item per ctp

"""

rule perm_age_labels:
    input: 
        sce_yng_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_yng-05",
        sce_old_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_old-05",
        lrdb = LRDB
    output:
        perm_scores_yng = OUTPUT_DAT + "/03_perm_age/perm_scores_{species}_yng",
        perm_scores_old = OUTPUT_DAT + "/03_perm_age/perm_scores_{species}_old"
    params:
        perm_functions = "../../source/cci_functions_perm.R",
        main_functions = "../../source/cci_functions_calculation.R",
        iterations = VALUES["iterations"],
        nr_cores = VALUES["nr_cores_age"],
        min_perc = VALUES["min_perc"],
        top_level = VALUES["top_level"],
        ages = ages
    script:
        "scripts/03_perm_age_labels.R"  
        
"""

Get permutation statistics

The actual scores from the original CCI objects are obtained and the test
statistic = score delta (score old - score yng) is calculated.

For each permuted set of of "yng" and "old" CCI scores, the score delta is 
calculated accordingly. 
This is calculated for each interaction and for each interacting cell type pair.
From the n permuted score delta sets, statistics are estimated or calculated.
Also, vectors of the n permuted score deltas per interaction per cell type pair
are exported for visualisation.

Output = list, one item per cell type pair. Each list contains 1 dataframe with 
statistics (res_df), and a df containing permuted score deltas (vis_df) for
visualisation.

"""

rule perm_age_stats:
    input: 
        cci_yng = OUTPUT_BASE + "/cci_objects/01_cci_preparation/09_intl/interaction_list_{species}_yng",
        cci_old = OUTPUT_BASE + "/cci_objects/01_cci_preparation/09_intl/interaction_list_{species}_old",
        perm_scores_yng = rules.perm_age_labels.output.perm_scores_yng,
        perm_scores_old = rules.perm_age_labels.output.perm_scores_old,
    output:
        res_lists = OUTPUT_DAT + "/04_pres_age/res_lists_ages_{species}"
    params:
        perm_functions = "../../source/cci_functions_perm.R",
        iterations = VALUES["iterations"],
        nr_cores = VALUES["nr_cores_age"],
        min_sp_age_perc = VALUES["min_sp_age_perc"]
    script:
        "scripts/04_perm_age_stats.R"
       
rule perm_age_report:
    input: 
        res_lists = rules.perm_age_stats.output,
        sce_yng_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_yng-05",
        sce_old_input = OUTPUT_BASE + "/cci_objects/01_cci_preparation/05_down/sce_{species}_old-05"
    output:
        OUTPUT_REP + "/report_perm_age/{species}_report.html"
    params:
        iterations = VALUES["iterations"]
    threads:
        4
    script:
        "perm_age_report.Rmd"  
        
rule perm_age_summary:
    input: 
        res_list_input_path = expand(OUTPUT_DAT + "/04_pres_age/res_lists_ages_{species}", species = species)
    output:
        summary_df = OUTPUT_DAT + "/05_psum_age/summary_df"
    params:
        species = species
    script:
        "scripts/05_perm_age_summary.R"

rule perm_age_summary_report:
    input: 
        summary_df = rules.perm_age_summary.output,
    output:
        OUTPUT_REP + "/report_perm_age/perm_summary.html"
    params:
        iterations = VALUES["iterations"],
        colors_path = COLORS,
        colors = "../../source/colors.R"
    threads:
        4
    script:
        "perm_age_summary.Rmd"  
