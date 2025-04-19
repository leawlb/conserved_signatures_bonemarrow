#!/bin/python 

import pandas as pd

#-------------------------------------------------------------------------------

OUTPUT_BASE = config["base"] + config["scRNAseq_data_paths"]["main"]

OUTPUT_DAT = OUTPUT_BASE + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own"
OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures/01_reclustering_own"

print(OUTPUT_DAT)
print(OUTPUT_REP)

COLORS = config["base"] + config["metadata_paths"]["colors"]

CELL_TYPES_EXCLUDE = config["values"]["03_sce_analysis"]["cell_types_exclude"]

# permutation: each gene set individually vs random 
RUN_SIGN_RAND_OWN_PERM = config["run_sign_rand_own_permutation"]
RUN_MARK_RAND_OWN_PERM = config["run_mark_rand_own_permutation"]
RUN_MMMS_RAND_OWN_PERM = config["run_mmms_rand_own_permutation"]

# permutation: marker sets vs signature + random
RUN_GNST_SIGN_OWN_PERM = config["run_genesets_sign_own_permutation"]

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"]
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"]
ENSEMBL_ZEB = config["base"] + config["metadata_paths"]["ensembl_zeb"]
ENSEMBL_NMR = config["base"] + config["metadata_paths"]["ensembl_nmr"]

#-------------------------------------------------------------------------------

METADATA = pd.read_csv(config["table"])
def get_list(metadata, column):
  values = METADATA[column]
  values = values.drop_duplicates()
  values = values.squeeze()
  values = values.tolist()
  return(values)

fractions = get_list(metadata = METADATA, column = "Fraction_ID")

#-------------------------------------------------------------------------------

targets = []
for f in fractions:
  targets = targets + [OUTPUT_DAT + "/01_gens/geneset_list_" + f]

  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_sign_" + f]
  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mark_" + f]
  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_ndge_" + f]
  targets = targets + [OUTPUT_DAT + "/02_endf/ensembl_mmms_" + f]
  
  targets = targets + [OUTPUT_DAT + "/03_recl/sce_" + f]
  targets = targets + [OUTPUT_DAT + "/04_rcls/score_df_" + f]
  targets = targets + [OUTPUT_REP + "/reclustering_own_report_" + f + ".html"]

  if RUN_SIGN_RAND_OWN_PERM:
    targets = targets + [OUTPUT_DAT + "/05_psig/perm_score_df_" + f]
    targets = targets + [OUTPUT_DAT + "/08_expp/sign-vs-rand_" + f]
    targets = targets + [OUTPUT_REP + "/perm_signature_" + f + ".html"]

  if RUN_MARK_RAND_OWN_PERM:
    targets = targets + [OUTPUT_DAT + "/05_pmrk/perm_score_df_" + f]
    targets = targets + [OUTPUT_DAT + "/08_expp/mark-vs-rand_" + f]

  if RUN_MMMS_RAND_OWN_PERM:
    targets = targets + [OUTPUT_DAT + "/05_pmmm/perm_score_df_" + f]
    targets = targets + [OUTPUT_DAT + "/08_expp/mmms-vs-rand_" + f]
    
    
  if RUN_GNST_SIGN_OWN_PERM:
    targets = targets + [OUTPUT_DAT + "/06_perg/perm_score_df_mark_" + f]
    targets = targets + [OUTPUT_DAT + "/06_perg/perm_score_df_mmms_" + f]
    targets = targets + [OUTPUT_DAT + "/06_perg/perm_score_df_mmms_mark_" + f]
    targets = targets + [OUTPUT_DAT + "/08_expp/mark-vs-signrand_" + f]
    targets = targets + [OUTPUT_DAT + "/08_expp/mmms-vs-signrand_" + f]
    #targets = targets + [OUTPUT_DAT + "/08_expp/mmms-vs-markrand_" + f] # not required atm, remove once decision is final
    targets = targets + [OUTPUT_REP + "/perm_genesets_" + f + ".html"]

#targets = targets + [OUTPUT_REP + "/genesets_summary.html"]

# downloaded code from mcclust, explanation below
targets = targets + ["../../source/mcclust/mcclust-master/R/vi.dist.R"]

#-------------------------------------------------------------------------------

localrules: all, download_mcclust_vidist_code, sign_vs_rand, mark_vs_rand, mmms_vs_rand, mark_vs_signrand, mmms_vs_signrand

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""   
# EXTRACT GENE SETS
# 
# Get different gene sets:
# - signature genes = conserved marker genes + nDGEs=stable genes (SIGN)
# - conserved marker genes (MARK)
# - all BL6 marker genes (MMMS)
# (- nDGEs=stable genes)
"""

# export list of gene sets
rule export_genesets:
    input:
        celltype_ndge_list = OUTPUT_BASE + "/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/06_nres/PC_0.05_FC_1.5/shared_genes_{fraction}_celltypes",
        marker_cons = OUTPUT_BASE + "/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_{fraction}.RData",
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
    output:
        geneset_list = OUTPUT_DAT + "/01_gens/geneset_list_{fraction}"
    resources:
          mem_mb=20000,
          queue = "medium-debian"
    threads: 4
    params:
        cts_exclude = CELL_TYPES_EXCLUDE
    script:
        "01_scripts_own/01_export_genesets.R"    

# visualise gene sets across species and cell types
rule genesets_summary:
    input: 
        geneset_list_hsc = OUTPUT_DAT + "/01_gens/geneset_list_hsc",
        geneset_list_str = OUTPUT_DAT + "/01_gens/geneset_list_str",
        sce_hsc = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_hsc-10",
        sce_str = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_str-10"
    output:
        OUTPUT_REP + "/genesets_summary.html"
    resources:
          mem_mb=35000,
          queue = "medium-debian"
    threads: 4
    params:
        cts_exclude = CELL_TYPES_EXCLUDE,
        colors_path = COLORS,
        functions = "../../source/sce_functions.R",
        functions_reclustering = "../../source/sce_functions_reclustering.R",
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "01_genesets_summary.Rmd"

# prepare ensembl conversion tables for each gene set to be tested
# ensembl datasets are downloaded in 00_prepare_datasets_snakefile.py
# includes species for other datasets in 02_signatures_reclustering_other
rule prepare_ensembl:
    input:
        geneset_list = rules.export_genesets.output,
        sce_inp = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    output:
        ensembl_sign = OUTPUT_DAT + "/02_endf/ensembl_sign_{fraction}",
        ensembl_mark = OUTPUT_DAT + "/02_endf/ensembl_mark_{fraction}",
        ensembl_ndge = OUTPUT_DAT + "/02_endf/ensembl_ndge_{fraction}",
        ensembl_mmms = OUTPUT_DAT + "/02_endf/ensembl_mmms_{fraction}"
    resources:
          mem_mb=20000,
          queue = "medium-debian"
    threads: 4
    script:
        "01_scripts_own/02_prepare_ensembl.R"
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
       
"""
# RECLUSTERING AND PERMUTATION TESTS
# 
# Reclustering datasets to test the ability of the different gene sets
# to capture cell identity.
# Using our own datasets with:
#
# - signature genes
# - conserved marker genes
# - all BL6 marker genes 
# - nDGEs=stable genes

# Important: subclustering genes are removed completely

# Additionally, permutation tests are performed using: 
# the same number of random genes for each gene set
# OR
# signature + required remaining amount of random genes for each marker set
"""

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# RE-CLUSTERING
#
# Data will be reclustered as the original dataset with some adjustments.
# - SCE object are subsetted to each gene set
# - PCA is performed again using non-batch corrected values

# Then, clustering is performed at identical resolution and k = number of 
# neighbors as before.
"""

# reclustering & calculating scores in total requires a specific function
# which is not easily available anymore anaconda distribution r channel)
# so in order to still be able to use it, all code from the corresponding
# github repository is downloaded into ../../source

# the function is called vi.dist and is from the mcclust package
# the entire github repository code will be downloaded and put into source/mcclust/
# https://doi.org/10.32614/CRAN.package.mcclust
# https://github.com/cran/mcclust
# Author: Arno Fritsch
# the exact download link is download_zip_url

# last downloaded 2025-01-22
# will not be downloaded again as long as mcclust-master.zip exists
rule download_mcclust_vidist_code:
    output:
        path_output = "../../source/mcclust/mcclust-master/R/vi.dist.R",
    params:
        download_zip_url = "https://github.com/cran/mcclust/archive/refs/heads/master.zip",
        output_directory = "../../source/mcclust" 
    conda:
        "../../envs/unzip.yml"
    shell:
        """
        mkdir -p {params.output_directory} 
        cd {params.output_directory} 
        if [ -e "mcclust-master.zip" ]; then
          echo 'File already exists' >&2
        else   
          curl {params.download_zip_url} -O -L -J 
        fi
          unzip -o mcclust-master.zip -d {params.output_directory}
        """


rule reclustering_own:
    input:
        sce_input = OUTPUT_BASE + "/sce_objects/02_sce_anno/10_anns/sce_{fraction}-10",
        geneset_list = rules.export_genesets.output
    params:
        k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
        resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
        nr_cores = 4, # 4 genesets total for reclustering
        cts_exclude = CELL_TYPES_EXCLUDE,
        functions_reclustering = "../../source/sce_functions_reclustering.R"
    output:
        sce_output = OUTPUT_DAT + "/03_recl/sce_{fraction}"
        # contains re-clustered vectors
        # cts_exclude already removed + factor levels removed
        # HSC sub-sampled to 25,000 cells
    resources:
          mem_mb=48000,
          queue = "medium-debian"
    threads: 4
    script:
        "01_scripts_own/03_reclustering_own.R"

# get the reclustering scores for our own re-clustered datasets
# the conda environment can be used to calculate many more types of scores
# than are currently required - but were only used for testing
rule reclustering_own_scores:
    input:
        sce_input = rules.reclustering_own.output,
    params:
        cts_exclude = CELL_TYPES_EXCLUDE,
        functions_reclustering = "../../source/sce_functions_reclustering.R"
    conda:
        "../../envs/reclust_scores_perm.yml"
    resources:
          mem_mb=8000,
          queue = "medium-debian"
    threads: 4
    output:
        score_df = OUTPUT_DAT + "/04_rcls/score_df_{fraction}"
    script:
        "01_scripts_own/04_reclustering_own_scores.R"
        
# visualise re-clustering of own datasets, including scores
rule reclustering_own_report:
    input: 
        sce_input = rules.reclustering_own.output,
        score_df = rules.reclustering_own_scores.output
    output:
        OUTPUT_REP + "/reclustering_own_report_{fraction}.html"
    resources:
          mem_mb=10000,
          queue = "medium-debian"
    threads: 4
    params:
        colors_path = COLORS,
        plotting = "../../source/plotting.R",
        colors = "../../source/colors.R"
    script:
        "01_reclustering_own_report.Rmd"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# PERMUTATION TESTS

# 1. Do the three gene sets  perform significantly better than random?
# Permutation tests to compare:
# - signature genes
# - conserved markers
# - BL6 markers (species-specific)
# to gene sets comprised of x random genes 

# 2. Do other gene sets perform better than signature genes?
#    Or is it just the larger number of genes?
# Permutation tests to compare:
# - conserved marker genes
# - all BL6 marker genes (species-specific)
# to gene sets comprised of signature genes + x random genes 
# or conserved markers genes + x random genes


# Re-cluster own datasets with random genes n = iteration times and compare
# the re-clustering scores of the original seurat object to the normal 
# distribution of the permuted scores.

# For permutation, the same resolution as the original gene set is chosen, 
# because each random set of genes might have a different optimal resolution. 

"""

#-------------------------------------------------------------------------------

# permutation: signature vs. random (same number)
if RUN_SIGN_RAND_OWN_PERM:

  rule perm_sign_rand_own:
      input:
          sce_input = rules.reclustering_own.output,
          # cts_exclude already excluded + factor levels removed
          # HSCs subsampled to 25,000 cells
          geneset_list = rules.export_genesets.output,
      params:
          k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          iterations = config["values"]["03_sce_analysis"]["iterations"],
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          functions_reclustering = "../../source/sce_functions_reclustering.R",
          cts_exclude = CELL_TYPES_EXCLUDE,       
          cons_level_use = "conserved_signature"      
      conda:
          "../../envs/reclust_scores_perm.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/05_psig/perm_score_df_{fraction}"
      resources:
          mem_mb=180000,
          queue = "long-debian"
      threads: config["values"]["03_sce_analysis"]["nr_cores"]
      script: 
          "01_scripts_own/05_permutation_rand.R"

  # visualise signature permutation
  rule perm_sign_rand_own_report:
      input:
          score_df = rules.reclustering_own_scores.output,
          perm_score_df = rules.perm_sign_rand_own.output
      params:
          cons_level_use = "conserved_signature"
      output:
          OUTPUT_REP + "/perm_signature_{fraction}.html"
      resources:
          mem_mb=1000,
          queue = "medium-debian"
      threads:4
      script: 
          "01_permutation_own_background_report.Rmd"

  # compare signature to random background
  rule sign_vs_rand:
      input:
          orig_score_df_input = rules.reclustering_own_scores.output.score_df,
          perm_score_df_input = rules.perm_sign_rand_own.output.perm_score_df
      params:
          cons_level_use = "conserved_signature",
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          comparison = "sign-vs-rand"
      output:
          pval_score_df_output = OUTPUT_DAT + "/08_expp/sign-vs-rand_{fraction}"
      #resources:
      #    mem_mb=50000  
      script:
          "01_scripts_own/08_export_pval.R"

if RUN_MARK_RAND_OWN_PERM:

  rule perm_mark_rand_own:
      input:
          sce_input = rules.reclustering_own.output,
          geneset_list = rules.export_genesets.output,
      params:
          k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          iterations = config["values"]["03_sce_analysis"]["iterations"],
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          functions_reclustering = "../../source/sce_functions_reclustering.R",
          cts_exclude = CELL_TYPES_EXCLUDE,       
          cons_level_use = "conserved_markers"      
      conda:
          "../../envs/reclust_scores_perm.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/05_pmrk/perm_score_df_{fraction}" ##### OUTPUT TEMP
      resources:
          mem_mb=180000,
          queue = "long-debian"
      threads: config["values"]["03_sce_analysis"]["nr_cores"] 
      script: 
          "01_scripts_own/05_permutation_rand.R"

  # compare conserved markers to random background
  rule mark_vs_rand:
      input:
          orig_score_df_input = rules.reclustering_own_scores.output.score_df,
          perm_score_df_input = rules.perm_mark_rand_own.output.perm_score_df
      params:
          cons_level_use = "conserved_markers",
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          comparison = "mark-vs-rand"
      output:
          pval_score_df_output = OUTPUT_DAT + "/08_expp/mark-vs-rand_{fraction}"
      #resources:
      #    mem_mb=50000  
      script:
          "01_scripts_own/08_export_pval.R"

if RUN_MMMS_RAND_OWN_PERM:

  rule perm_mmms_rand_own:
      input:
          sce_input = rules.reclustering_own.output,
          geneset_list = rules.export_genesets.output,
      params:
          k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"],
          iterations = config["values"]["03_sce_analysis"]["iterations"],
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          functions_reclustering = "../../source/sce_functions_reclustering.R",
          cts_exclude = CELL_TYPES_EXCLUDE,       
          cons_level_use = "mmusall_markers"      
      conda:
          "../../envs/reclust_scores_perm.yml"
      output:
          perm_score_df = OUTPUT_DAT + "/05_pmmm/perm_score_df_{fraction}" ##### OUTPUT TEMP
      resources:
          mem_mb=180000,
          queue = "long-debian"
      threads: config["values"]["03_sce_analysis"]["nr_cores"]  
      script: 
          "01_scripts_own/05_permutation_rand.R"

  # compare BL6 markers to random background
  rule mmms_vs_rand:
      input:
          orig_score_df_input = rules.reclustering_own_scores.output.score_df,
          perm_score_df_input = rules.perm_mmms_rand_own.output.perm_score_df
      params:
          cons_level_use = "mmusall_markers",
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          comparison = "mmms-vs-rand"
      output:
          pval_score_df_output = OUTPUT_DAT + "/08_expp/mmms-vs-rand_{fraction}"
      #resources:
      #    mem_mb=50000  
      script:
          "01_scripts_own/08_export_pval.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# permutation: marker genes vs. signature + random (same number)
if RUN_GNST_SIGN_OWN_PERM:

  rule perm_genesets_signature_rand_own:
      input:
          sce_input = rules.reclustering_own.output,
          # cts_exclude already excluded + factor levels removed
          # HSCs subsampled to 25,000 cells
          geneset_list = rules.export_genesets.output,
      params:
          functions_reclustering = "../../source/sce_functions_reclustering.R",
          k_graph_list = config["values"]["02_sce_anno"]["k_graph_list"],
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          nr_cores = config["values"]["03_sce_analysis"]["nr_cores"],
          iterations = config["values"]["03_sce_analysis"]["iterations"],
          cut_off_prop = config["values"]["03_sce_analysis"]["reclustering_cutoff_prop"]
      conda:
          "../../envs/reclust_scores_perm.yml"
      output:
          perm_score_df_mark = OUTPUT_DAT + "/06_perg/perm_score_df_mark_{fraction}",
          perm_score_df_mmms = OUTPUT_DAT + "/06_perg/perm_score_df_mmms_{fraction}",
          perm_score_df_mmms_mark = OUTPUT_DAT + "/06_perg/perm_score_df_mmms_mark_{fraction}"
      resources:
          mem_mb=180000,
          queue = "long-debian"
      threads: config["values"]["03_sce_analysis"]["nr_cores"]  
      script: 
          "01_scripts_own/06_permutation_genesets_signrand.R"

  # visualise reclustering permutation
  rule perm_genesets_signature_rand_own_report:
      input:
          score_df = rules.reclustering_own_scores.output,
          geneset_list = rules.export_genesets.output,
          perm_score_df_mark = rules.perm_genesets_signature_rand_own.output.perm_score_df_mark,
          perm_score_df_mmms = rules.perm_genesets_signature_rand_own.output.perm_score_df_mmms,
          perm_score_df_mmms_mark = rules.perm_genesets_signature_rand_own.output.perm_score_df_mmms_mark
      output:
          OUTPUT_REP + "/perm_genesets_{fraction}.html"
      resources:
          mem_mb=1000,
          queue = "medium-debian"
      threads: 4
      script: 
          "01_permutation_own_genesets_report.Rmd"
    
  # compare conserved markers to signature + random
  rule mark_vs_signrand:
      input:
          orig_score_df_input = rules.reclustering_own_scores.output.score_df,
          perm_score_df_input = rules.perm_genesets_signature_rand_own.output.perm_score_df_mark
      params:
          cons_level_use = "conserved_markers",
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          comparison = "mark-vs-signrand"
      output:
          pval_score_df_output = OUTPUT_DAT + "/08_expp/mark-vs-signrand_{fraction}"
      #resources:
      #    mem_mb=50000  
      script:
          "01_scripts_own/08_export_pval.R"
          
  # compare mmusall_markers (all BL6 markers) to signature + random
  rule mmms_vs_signrand:
      input:
          orig_score_df_input = rules.reclustering_own_scores.output.score_df,
          perm_score_df_input = rules.perm_genesets_signature_rand_own.output.perm_score_df_mmms
      params:
          cons_level_use = "mmusall_markers",
          resolution_louvain_list = config["values"]["02_sce_anno"]["resolution_louvain_list"],
          comparison = "mmms-vs-signrand"
      output:
          pval_score_df_output = OUTPUT_DAT + "/08_expp/mmms-vs-signrand_{fraction}"
      #resources:
      #    mem_mb=50000  
      script:
          "01_scripts_own/08_export_pval.R"
          
  """
  # At the moment, this step is not necessary.  
  # Remove from 01_permutation_own_genesets_report once decision is final.
  # compare mmusall_markers (all BL6 markers) to conserved markers + random
  rule mmms_vs_markrand:
      input:
          orig_score_df_input = rules.reclustering_own_scores.output.score_df,
          perm_score_df_input = rules.permutation_genesets.output.perm_score_df_mmms_mark
      params:
          cons_level_use = "mmusall_markers",
          comparison = "mmms-vs-markrand"
      output:
          pval_score_df_output = OUTPUT_DAT + "/08_expp/mmms-vs-markrand_{fraction}"
      script:
          "01_scripts_own/08_export_pval.R"
  """
  
