#!/bin/python 

# import pandas as pd

# the purpose of this snakefile is to
# - prepare publically available scRNAseq datasets and for re-clustering
# - download ensembl ID conversion tables between MMUS and other species

# all of the generated data is saved in METADATA

# last run 2024-08-07

#-------------------------------------------------------------------------------

DIR_RECLUSTERING = config["base"] + config["metadata_paths"]["reclustering"]

OUTPUT_REP = OUTPUT_BASE + "/sce_objects/reports/03_sce_analysis/04_signatures"

RAN_LI_IPYNB = config["ran_li_ipynb"]

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"]
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"]
ENSEMBL_ZEB = config["base"] + config["metadata_paths"]["ensembl_zeb"]
ENSEMBL_NMR = config["base"] + config["metadata_paths"]["ensembl_nmr"]

#-------------------------------------------------------------------------------

targets = []

# ensembl conversion sheets
targets = targets + [ENSEMBL_MUS]
targets = targets + [ENSEMBL_HUM]
targets = targets + [ENSEMBL_ZEB]
targets = targets + [ENSEMBL_NMR]

# downloaded datasets (tabula sapiens)
targets = targets + [DIR_RECLUSTERING + "/raw_hum/ts_bone_marrow"]
targets = targets + [DIR_RECLUSTERING + "/raw_hum/ts_hscs_progenitors"]
targets = targets + [DIR_RECLUSTERING + "/raw_hum/ts_all_stromal"]

# prepared datasets (tabula sapiens)
targets = targets + [DIR_RECLUSTERING + "/prepared/ts_bone_marrow"]
targets = targets + [DIR_RECLUSTERING + "/prepared/ts_hscs_progenitors"]
targets = targets + [DIR_RECLUSTERING + "/prepared/ts_all_stromal"]

# prepared dataset (Li) only after ipynb script was executed for pre-processing
if RAN_LI_IPYNB:
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_only_counts"]
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_nocounts"]
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_all"]
  targets = targets + [DIR_RECLUSTERING + "/Li/sce_li_stromal"]
  targets = targets + [DIR_RECLUSTERING + "/prepared/li_all_stromal"]

# download exotic datasets
targets = targets + [DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_sorted_BM.RData"]
targets = targets + [DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_whole_BM.RData"]
targets = targets + [DIR_RECLUSTERING + "/raw_zeb/zeb_clusters.tsv"]

# download and unzip zebrafish datasets
# upon unzipping, files retain the names of the assays given by authors, 
# so I am using to the given file names as wildcards
endings = [".mtx", ".mtx_cols", ".mtx_rows"]
# this is what the unzipped files are called
assays = ["aggregated_filtered_normalised_counts", "aggregated_filtered_counts", "expression_tpm"]

for a in assays:
  targets = targets + [DIR_RECLUSTERING + "/raw_zeb/" + a + ".zip"]
  for e in endings:
    targets = targets + [DIR_RECLUSTERING + "/raw_zeb/" + a + "/E-MTAB-5530." + a + e]

# prepare exotic files
targets = targets + [DIR_RECLUSTERING + "/prepared/zeb_all_hspc"]
targets = targets + [DIR_RECLUSTERING + "/prepared/nmr_sorted_hspc"]
targets = targets + [DIR_RECLUSTERING + "/prepared/nmr_whole_hspc"]

#-------------------------------------------------------------------------------

wildcard_constraints: 
    fraction="[a-z]+"

localrules: all, unzip_zebrafish

rule all: 
    input:
        targets
        
#-------------------------------------------------------------------------------

# download ensembl conversion tables 
rule download_ensembl:
    output:
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    script:
        "prepare_datasets/download_ensembl.R"

#-------------------------------------------------------------------------------
"""
Human
"""

# download datasets
# the directories "/Li/download" and "/Li/data" must be generated manually
rule download_ts_datasets:
    output:
        tabula_sapiens_bone_marrow = DIR_RECLUSTERING + "/raw_hum/ts_bone_marrow",
        tabula_sapiens_hsc_progenitors = DIR_RECLUSTERING + "/raw_hum/ts_hscs_progenitors",
        tabula_sapiens_stromal = DIR_RECLUSTERING + "/raw_hum/ts_all_stromal"
    params:
        li_bone_marrow = DIR_RECLUSTERING + "/Li"
    script:
        "prepare_datasets/download_ts_files.R"
        
# prepare datasets - Tabula Sapiens
rule prepare_ts_datasets:
    input:
        ts_bone_marrow_input = rules.download_ts_datasets.output.tabula_sapiens_bone_marrow,
        ts_hsc_progenitors_input = rules.download_ts_datasets.output.tabula_sapiens_hsc_progenitors,
        ts_stromal_input = rules.download_ts_datasets.output.tabula_sapiens_stromal
    output:
        ts_bone_marrow_output = DIR_RECLUSTERING + "/prepared/ts_bone_marrow",
        ts_hsc_progenitors_output = DIR_RECLUSTERING + "/prepared/ts_hscs_progenitors",
        ts_stromal_output = DIR_RECLUSTERING + "/prepared/ts_all_stromal"
    script:
        "prepare_datasets/prepare_ts_datasets.R"

"""
# Manually pre-process Li et al. human stromal bone marrow dataset
# https://doi.org/10.7554/elife.81656
# according to https://github.com/Hongzhe2022/MSC_BM_scripts.

# I also added this chunk at the end of the .ipynb script so that adata with 
# clustering information is saved:

  # adding one bit to save with the publication cluster names:
  # this overwrites the previously saved version of adata

  if os.path.exists(ofile):
      os.remove(ofile)
  adata.write(ofile)
  print (ofile)

# After manual pre-processing, convert to Seurat and prepare for analysis.
# Only execute this rule when Load_GEO_data_analyze_it.ipynb has been executed.
"""
if RAN_LI_IPYNB:
  rule prepare_li:
      input:
          h5ad_input = DIR_RECLUSTERING + "/Li/download/GSE190965_analyzed.h5ad",
          h5ad_raw_input = DIR_RECLUSTERING + "/Li/download/GSE190965_raw.h5ad"
      output:
          sce_raw_output_1 = DIR_RECLUSTERING + "/Li/sce_li_only_counts",
          sce_output_1 = DIR_RECLUSTERING + "/Li/sce_li_nocounts",
          sce_output_2 = DIR_RECLUSTERING + "/Li/sce_li_all",
          sce_output_3 = DIR_RECLUSTERING + "/Li/sce_li_stromal",
          seurat_output = DIR_RECLUSTERING + "/prepared/li_all_stromal"
      conda:
          "../../envs/zellkonverter_li.yml"
      script:
          "prepare_datasets/convert_lidata.R"
        
#-------------------------------------------------------------------------------
"""
Mouse
"""

# download harris atlas loom file
rule download_harris_atlas:
    input:
        "https://labshare.cshl.edu/shares/gillislab/resource/HSC_atlas/processed_droplet_data.loom"
    output:
        DIR_RECLUSTERING + "/raw_mus/harris_atlas/processed_droplet_data.loom"
    shell:
        """
        wget {input} 
        """

# TABULA MURIS
# download the zipped folder with counts csv
# other data are downloaded in download_tm_R
rule download_tm_zipped:
    input:
        "https://figshare.com/ndownloader/files/10700143"
    output:
        DIR_RECLUSTERING + "/raw_mus/tabula_muris/FACS/Marrow-counts.csv"
    shell:
        """
        wget {input} 
        """

#-------------------------------------------------------------------------------
"""
NMR & Zebrafish
"""

rule download_exotic_datasets:
    output:
        Seurat_hgl_sorted_BM = DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_sorted_BM.RData",
        Seurat_hgl_whole_BM = DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_whole_BM.RData",
        zeb_clustering_file = DIR_RECLUSTERING + "/raw_zeb/zeb_clusters.tsv",
        expression_tpm = DIR_RECLUSTERING + "/raw_zeb/expression_tpm.zip",
        aggregated_filtered_normalised_counts = DIR_RECLUSTERING + "/raw_zeb/aggregated_filtered_normalised_counts.zip",
        aggregated_filtered_counts = DIR_RECLUSTERING + "/raw_zeb/aggregated_filtered_counts.zip"
    script:
        "prepare_datasets/download_exotic_datasets.R"  

rule unzip_zebrafish:
    input:
        DIR_RECLUSTERING + "/raw_zeb/{assay}.zip",
    output:
        expand(DIR_RECLUSTERING + "/raw_zeb/{{assay}}/E-MTAB-5530.{{assay}}{ending}", ending = endings)
    params:
        output_directory = DIR_RECLUSTERING + "/raw_zeb/{assay}/"
    conda:
        "../../envs/unzip.yml"
    shell:
        """
        mkdir -p {params.output_directory} 
        unzip {input} -d {params.output_directory}
        """

rule prepare_exotic_datasets:
    input:
        assays = expand(rules.unzip_zebrafish.output, assay = assays),
        zeb_clustering_file = rules.download_exotic_datasets.output.zeb_clustering_file,
        Seurat_hgl_sorted_BM = rules.download_exotic_datasets.output.Seurat_hgl_sorted_BM,
        Seurat_hgl_whole_BM = rules.download_exotic_datasets.output.Seurat_hgl_whole_BM
    output:
        zeb_seu = DIR_RECLUSTERING + "/prepared/zeb_all_hspc",
        nmr_sorted_hspc = DIR_RECLUSTERING + "/prepared/nmr_sorted_hspc",
        nmr_whole_hspc = DIR_RECLUSTERING + "/prepared/nmr_whole_hspc"
    script:
        "prepare_datasets/prepare_exotic_datasets.R"
        
#-------------------------------------------------------------------------------

rule report_datasets:
    input:
        zeb_seu = rules.prepare_exotic_datasets.output.zeb_seu,
        nmr_sorted_hspc = rules.prepare_exotic_datasets.output.nmr_sorted_hspc,
        li_all_stromal = rule.prepare_li.output.seurat_output,
        ts_bone_marrow = rules.prepare_ts_datasets.output,ts_bone_marrow_output,
        ts_hsc_progenitors = rules.prepare_ts_datasets.output,ts_hsc_progenitors_output,
        ts_stromal = rules.prepare_ts_datasets.output,ts_stromal_output
    output:
        OUTPUT_REP + "/report_other_datasets.html"
    script:
        "report_other_datasets.Rmd"

