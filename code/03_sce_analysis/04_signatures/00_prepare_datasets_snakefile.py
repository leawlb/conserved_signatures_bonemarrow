#!/bin/python 

# last run completely 2025-04-17
# including downloads on that same date (except for NMR and zebrafish)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""

The purpose of this snakefile is to:
 - download ensembl ID conversion tables between MMUS and other species
 - prepare publically available scRNAseq datasets for re-clustering

All of the data generated from this snakefile will be saved as METADATA in 
DIR_RECLUSTERING (see below).

NMR has been excluded from this analysis.




How to use this snakefile?

1. Upon first run, manually generate the directories:

  BASE + "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/"
  BASE + "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/Li/download"
  BASE + "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/Li/data"
  BASE + "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/FACS/"
  BASE + "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/weinreb/"


2. Some data must be downloaded and prepared manually (Li et al.):

  Publication: https://doi.org/10.7554/eLife.81656
  Github: https://github.com/Hongzhe2022/MSC_BM_scripts

  From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190965
    GSE190965_barcodes.tsv.gz (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Fbarcodes%2Etsv%2Egz)
    GSE190965_features.tsv.gz (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Ffeatures%2Etsv%2Egz)
    GSE190965_spliced.mtx.gz (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Fspliced%2Emtx%2Egz)
    GSE190965_unspliced.mtx.gz (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Funspliced%2Emtx%2Egz)

  --> put into "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/Li/download"


  Download the cluster information from their github at:
  https://github.com/Hongzhe2022/MSC_BM_scripts/blob/master/data/Hongzhe_clustering_and_2dUMAP_louvain11_merged_95_renamed.csv
  
  --> put into "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/Li/data"


  Download the jupyter notebook from their github at:
  https://github.com/Hongzhe2022/MSC_BM_scripts/blob/master/scripts/Load_GEO_data_analyze_it.ipynb

  --> put into "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/Li/download"


3. Manually preprocess data (Li et al.) with their provided jupyter notebook

  Add this chunk of code to the end of the notebook to save the publication
  cluster names:
  
if os.path.exists(ofile):
  os.remove(ofile)
adata.write(ofile)
print (ofile)

  Run the entire jupyter notebook in "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/Li/download".
  
  For the appropriate kernel, install environment load_analyze_export_GEO.yml in /code/envs, activate it, 
  then install as kernel using:
  python -m ipykernel install --user --name=load_analyze_export_GEO

  There is a DeprecationWarning starting with scvelo 0.4.0 but the scvelo version in the given .yml file is 0.3.3 so it should be fine.
  There is a DeprecationWarning: is_categorical_dtype is deprecated and will be removed in a future version. (deprecated since pandas version 2.2.0)  but the pandas version in the given .yml file is 2.1.4 so it should be fine.
  
  --> The output (processed adata file) is "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/Li/download/GSE190965_analyzed.h5ad"


4. Some more data must be downloaded manually (Tabula muris)

  Publication: https://doi.org/10.1038/s41586-018-0590-4

  Download FACS.zip manually from:
  https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687

  --> put all contents into "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/FACS/" and unzip


5. Some more data must be downloaded manually (Tikhonova et al.)

  This file should already be downloaded and prepared for reference annotation
  see file /data/metadata/scRNAseq/01_sce_prep/references_raw/prep_annotation_references
  
  This object contains data from three prior publications, of which
  the merged object was used for rough reference annotation, and two subsets 
  from prior publicatons will be separately reclustered.
  
  "/data/metadata/scRNAseq/01_sce_prep/references_raw/bone-marrow-seurat.rds"

6. Also get mcclust (https://github.com/cran/mcclust) and put in code/source/

  see 01_reclustering_own_snakefile.py, automatically downloaded from there
  
"""


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# define paths

DIR_RECLUSTERING = config["base"] + config["metadata_paths"]["reclustering"]
OUTPUT_REP = config["base"] + config["scRNAseq_data_paths"]["main"] + "/sce_objects/reports/03_sce_analysis/04_signatures"

RAN_LI_IPYNB = config["ran_li_ipynb"]

ENSEMBL_MUS = config["base"] + config["metadata_paths"]["ensembl_mus"] 
ENSEMBL_HUM = config["base"] + config["metadata_paths"]["ensembl_hum"] 
# although not currently in use, download and add to ensembl table as well
# for the sake of keeping the original ensembl table
ENSEMBL_ZEB = config["base"] + config["metadata_paths"]["ensembl_zeb"] 
ENSEMBL_NMR = config["base"] + config["metadata_paths"]["ensembl_nmr"] 


print(DIR_RECLUSTERING)
print(OUTPUT_REP)
print(ENSEMBL_MUS)
print(ENSEMBL_HUM)
print(ENSEMBL_ZEB)
print(ENSEMBL_NMR)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# define targets

targets = []

#-------------------------------------------------------------------------------


# ensembl conversion sheets
targets = targets + [ENSEMBL_MUS]
targets = targets + [ENSEMBL_HUM]
targets = targets + [ENSEMBL_ZEB]
targets = targets + [ENSEMBL_NMR]

#-------------------------------------------------------------------------------

# human

# downloaded human datasets (tabula sapiens)
targets = targets + [DIR_RECLUSTERING + "/raw_hum/ts_bone_marrow"]
targets = targets + [DIR_RECLUSTERING + "/raw_hum/ts_hscs_progenitors"]
targets = targets + [DIR_RECLUSTERING + "/raw_hum/ts_all_stromal"]

# prepared datasets human (tabula sapiens)
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

#-------------------------------------------------------------------------------

# mouse

# downloaded mouse datasets (tabula muris and weinreb)
targets = targets + [DIR_RECLUSTERING + "/raw_mus/tabula_muris/metadata_FACS.csv"]
targets = targets + [DIR_RECLUSTERING + "/raw_mus/tabula_muris/annotations_facs.csv"]
targets = targets + [DIR_RECLUSTERING + "/raw_mus/weinreb/GSE140802_RAW.tar"]

# prepared mouse datasets (tabula muris and weinreb)

weinreb_outputs = [
  "GSM4185642_stateFate_inVitro_cell_barcodes.txt.gz",
  "GSM4185642_stateFate_inVitro_clone_matrix.mtx.gz",
  "GSM4185642_stateFate_inVitro_gene_names.txt.gz",
  "GSM4185642_stateFate_inVitro_library_names.txt.gz",
  "GSM4185642_stateFate_inVitro_metadata.txt.gz",
  "GSM4185642_stateFate_inVitro_normed_counts.mtx.gz",
  "GSM4185643_stateFate_inVivo_cell_barcodes.txt.gz",
  "GSM4185643_stateFate_inVivo_clone_matrix.mtx.gz",
  "GSM4185643_stateFate_inVivo_gene_names.txt.gz",
  "GSM4185643_stateFate_inVivo_library_names.txt.gz",
  "GSM4185643_stateFate_inVivo_metadata.txt.gz",
  "GSM4185643_stateFate_inVivo_normed_counts.mtx.gz",
  "GSM4185644_stateFate_cytokinePerturbation_cell_barcodes.txt.gz",
  "GSM4185644_stateFate_cytokinePerturbation_clone_matrix.mtx.gz",
  "GSM4185644_stateFate_cytokinePerturbation_gene_names.txt.gz",
  "GSM4185644_stateFate_cytokinePerturbation_library_names.txt.gz",
  "GSM4185644_stateFate_cytokinePerturbation_metadata.txt.gz",
  "GSM4185644_stateFate_cytokinePerturbation_normed_counts.mtx.gz",
  "GSM4185645_GMP_MDP_culture_cell_barcodes.txt.gz",
  "GSM4185645_GMP_MDP_culture_gene_names.txt.gz",
  "GSM4185645_GMP_MDP_culture_library_names.txt.gz",
  "GSM4185645_GMP_MDP_culture_normed_counts.mtx.gz"
]

for o in weinreb_outputs:
  targets = targets + [DIR_RECLUSTERING + "/raw_mus/weinreb/" + o]

targets = targets + [DIR_RECLUSTERING + "/prepared/mus_tm_bonemarrow"]
targets = targets + [DIR_RECLUSTERING + "/prepared/mus_weinreb_hspc"]
 
targets = targets + [DIR_RECLUSTERING + "/prepared/mus_tik_stromal"]
targets = targets + [DIR_RECLUSTERING + "/prepared/mus_bar_stromal"]


#-------------------------------------------------------------------------------

# exotic species
# zebrafish is currently needed
# nmr is not --> but they are simply downloaded together

# zebrafish files retain the names of the assays given by authors upon unzipping
# so need to use given file names as wildcards
endings = [".mtx", ".mtx_cols", ".mtx_rows"]
# this is what the unzipped files are called
assays = ["aggregated_filtered_normalised_counts", "aggregated_filtered_counts", "expression_tpm"]

# downloaded exotic datasets (zebrafish and naked mole rat)
targets = targets + [DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_sorted_BM.RData"]
targets = targets + [DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_whole_BM.RData"]
targets = targets + [DIR_RECLUSTERING + "/raw_zeb/zeb_clusters.tsv"]

for a in assays:
  targets = targets + [DIR_RECLUSTERING + "/raw_zeb/" + a + ".zip"]
  for e in endings:
    targets = targets + [DIR_RECLUSTERING + "/raw_zeb/" + a + "/E-MTAB-5530." + a + e]

# prepared exotic files
targets = targets + [DIR_RECLUSTERING + "/prepared/zeb_all_hspc"]
targets = targets + [DIR_RECLUSTERING + "/prepared/nmr_sorted_hspc"]
targets = targets + [DIR_RECLUSTERING + "/prepared/nmr_whole_hspc"]


# report
if RAN_LI_IPYNB:
  targets = targets + [OUTPUT_REP + "/prepare_datasets_report.html"]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

wildcard_constraints: 
    fraction="[a-z]+"

localrules: all, untar_weinreb, download_ensembl#download_ts_datasets, download_mouse_datasets, download_exotic_datasets, unzip_zebrafish,

rule all: 
    input:
        targets
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# download ensembl conversion tables for all species
# zebrafish and naked mole rat are not currently in use, but also downloaded
# and added to ensembl table for the sake of keeping up original processing

rule download_ensembl:
    output:
        ensembl_mus = ENSEMBL_MUS,
        ensembl_hum = ENSEMBL_HUM,
        ensembl_zeb = ENSEMBL_ZEB,
        ensembl_nmr = ENSEMBL_NMR
    resources:
        mem_mb=2000
    script:
        "00_prepare_datasets/download_ensembl.R"

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


"""
# Human
"""

# download tabula sapiens datasets
# download links are in the script
rule download_ts_datasets:
    output:
        tabula_sapiens_bone_marrow = DIR_RECLUSTERING + "/raw_hum/ts_bone_marrow",
        tabula_sapiens_hsc_progenitors = DIR_RECLUSTERING + "/raw_hum/ts_hscs_progenitors",
        tabula_sapiens_stromal = DIR_RECLUSTERING + "/raw_hum/ts_all_stromal"
    resources:
        mem_mb=500
    params:
        li_bone_marrow = DIR_RECLUSTERING + "/Li"
    script:
        "00_prepare_datasets/download_ts_files.R"
        
# prepare datasets - Tabula Sapiens
rule prepare_ts_datasets:
    input:
        ts_bone_marrow_input = rules.download_ts_datasets.output.tabula_sapiens_bone_marrow,
        ts_hsc_progenitors_input = rules.download_ts_datasets.output.tabula_sapiens_hsc_progenitors,
        ts_stromal_input = rules.download_ts_datasets.output.tabula_sapiens_stromal
    output:
        ts_bone_marrow = DIR_RECLUSTERING + "/prepared/ts_bone_marrow",
        ts_hscs_progenitors = DIR_RECLUSTERING + "/prepared/ts_hscs_progenitors",
        ts_all_stromal = DIR_RECLUSTERING + "/prepared/ts_all_stromal"
    resources:
        mem_mb=20000
    script:
        "00_prepare_datasets/prepare_ts_datasets.R"


"""
                                  # !!!!! # 
                                  
# The directories "/Li/download" and "/Li/data" must be generated manually!
# Files must be downloaded manually as descrived above! 
# Data must be pre-processed using their jupyter notebook!

                                  # !!!!! # 

# After manual pre-processing, convert analyzed adata object to Seurat and 
# prepare for analysis.

# Only execute this rule when Load_GEO_data_analyze_it.ipynb has been executed!!

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
          li_all_stromal = DIR_RECLUSTERING + "/prepared/li_all_stromal"
      conda:
          "../../envs/zellkonverter_li.yml"
      resources:
        mem_mb=40000
      script:
          "00_prepare_datasets/convert_lidata.R"
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# Mouse
"""

rule download_mouse_datasets:
    output:
        path_metadata_all = DIR_RECLUSTERING + "/raw_mus/tabula_muris/metadata_FACS.csv",
        path_annotation_all = DIR_RECLUSTERING + "/raw_mus/tabula_muris/annotations_facs.csv",
        path_weinreb = DIR_RECLUSTERING + "/raw_mus/weinreb/GSE140802_RAW.tar"
    resources:
        mem_mb=500
    script:
        "00_prepare_datasets/download_mouse_datasets.R"

# extract all weinreb download folders from .tar, but only need GSE4185642 files
rule untar_weinreb:
    input:
        rules.download_mouse_datasets.output.path_weinreb,
    output:
        expand(DIR_RECLUSTERING + "/raw_mus/weinreb/{output}", output = weinreb_outputs)
    params:
        output_folder = DIR_RECLUSTERING + "/raw_mus/weinreb/"
    shell:
        """
        cd {params.output_folder}
        tar -xvf {input}
        """


"""
                                  # !!!!! # 

# Requires manual download of the zipped "FACS" folder for tabula muris!

                                  # !!!!! # 
"""  
      
# prepare tabula muris dataset and weinreb dataset
rule prepare_mouse_datasets_hspc:
    input:
        path_counts_marrow = DIR_RECLUSTERING + "/raw_mus/tabula_muris/FACS/Marrow-counts.csv", # manually downloaded
        path_metadata_all = rules.download_mouse_datasets.output.path_metadata_all,
        path_annotation_all = rules.download_mouse_datasets.output.path_annotation_all,
        path_weinreb_base = DIR_RECLUSTERING + "/raw_mus/weinreb/"
    output:
        mus_tm_bonemarrow = DIR_RECLUSTERING + "/prepared/mus_tm_bonemarrow",
        mus_weinreb_hspc = DIR_RECLUSTERING + "/prepared/mus_weinreb_hspc"
    resources:
        mem_mb=15000
    script:
        "00_prepare_datasets/prepare_mouse_datasets_hspc.R"
        
# prepare tikhonova and baryawno dataset
# use the merged dataset from dolgalev reference dataset as input
rule prepare_mouse_datasets_stromal:
    input:
        seu_mus_str_ref = config["base"] + "/data/metadata/scRNAseq/01_sce_prep/references_raw/bone-marrow-seurat.rds" # manually downloaded
    output:
        mus_tik_stromal = DIR_RECLUSTERING + "/prepared/mus_tik_stromal",
        mus_bar_stromal = DIR_RECLUSTERING + "/prepared/mus_bar_stromal"
    resources:
        mem_mb=4000
    script:
        "00_prepare_datasets/prepare_mouse_datasets_stromal.R"
        
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# Exotic species

# Zebrafish is needed
# NMR is not needed

# both are still downloaded together

"""

rule download_exotic_datasets:
    output:
        Seurat_hgl_sorted_BM = DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_sorted_BM.RData",
        Seurat_hgl_whole_BM = DIR_RECLUSTERING + "/raw_nmr/Seurat_hgl_whole_BM.RData",
        zeb_clustering_file = DIR_RECLUSTERING + "/raw_zeb/zeb_clusters.tsv",
        expression_tpm = DIR_RECLUSTERING + "/raw_zeb/expression_tpm.zip",
        aggregated_filtered_normalised_counts = DIR_RECLUSTERING + "/raw_zeb/aggregated_filtered_normalised_counts.zip",
        aggregated_filtered_counts = DIR_RECLUSTERING + "/raw_zeb/aggregated_filtered_counts.zip"
    resources:
        mem_mb=500
    script:
        "00_prepare_datasets/download_exotic_datasets.R"  

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
        zeb_all_hspc = DIR_RECLUSTERING + "/prepared/zeb_all_hspc",
        nmr_sorted_hspc = DIR_RECLUSTERING + "/prepared/nmr_sorted_hspc",
        nmr_whole_hspc = DIR_RECLUSTERING + "/prepared/nmr_whole_hspc"
    resources:
        mem_mb=10000
    script:
        "00_prepare_datasets/prepare_exotic_datasets.R"


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

"""
# Report
"""

if RAN_LI_IPYNB:
  rule report_datasets:
      input:
          li_all_stromal = rules.prepare_li.output.li_all_stromal,
          ts_all_stromal = rules.prepare_ts_datasets.output.ts_all_stromal,
          ts_bone_marrow = rules.prepare_ts_datasets.output.ts_bone_marrow,
          ts_hscs_progenitors = rules.prepare_ts_datasets.output.ts_hscs_progenitors,
          mus_tm_bonemarrow = rules.prepare_mouse_datasets_hspc.output.mus_tm_bonemarrow,
          mus_weinreb_hspc = rules.prepare_mouse_datasets_hspc.output.mus_weinreb_hspc,
          mus_tik_stromal = rules.prepare_mouse_datasets_stromal.output.mus_tik_stromal,
          mus_bar_stromal = rules.prepare_mouse_datasets_stromal.output.mus_bar_stromal,
          zeb_all_hspc = rules.prepare_exotic_datasets.output.zeb_all_hspc,
          #nmr_sorted_hspc = rules.prepare_exotic_datasets.output.nmr_sorted_hspc
      params:
          reclustering_functions = "../../source/sce_functions_reclustering.R",
      resources:
          mem_mb=20000
      output:
          OUTPUT_REP + "/prepare_datasets_report.html"
      script:
          "00_prepare_datasets_report.Rmd"
