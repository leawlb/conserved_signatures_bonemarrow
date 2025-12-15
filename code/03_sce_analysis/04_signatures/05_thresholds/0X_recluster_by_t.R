
# recluster datasets from other species using different gene sets
# obtained from sliding the pvalue threshold of the FindMarkers analysis.
# use the same resolution as determined for the signature genes.

# !!!!
# this script is CUSTOM-WRITTEN for ts_hsc_progenitors, so it will not work
# for other datasets without manual adjustment of the resolution, ensembl IDs 
# and so on
# !!!!

# determine random number generator for sample
library(parallel)
RNGkind("L'Ecuyer-CMRG") # using this for parallel is necessary

library(tidyverse)
library(Seurat, quietly = TRUE)
library(SeuratObject, quietly = TRUE)

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# load objects from different thresholding

threshold_path <- snakemake@params[["threshold_path"]]
#threshold_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/05_thresholds/"
print(threshold_path)

# list of signature genes obtained when sliding pval threshold for nDGE analysis
signature_genes_by_t <- readRDS(base::paste0(
  threshold_path,
  "/signature_genes_by_t.rds"))

# lists of unique HSPC signature genes when sliding p value threshold 
# on nDEG analysis
names(signature_genes_by_t)
lapply(signature_genes_by_t, function(vec){print(length(vec))})

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# load other objects 

# I select ts_hspcs_progenitors, because it looks relatively stable in
# figure S4 even with different numbers of genes (no big jumps or bumps
# of scores with changes in resolution, because it is one of the first 
# human datasets and because the annotation has a better resolution than 
# the ts adult whole bone marrow

seu_dataset <- readRDS(snakemake@input[["seu_input"]])
#seu_dataset <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/prepared/ts_hscs_progenitors")

dataset_curr <- snakemake@params[["dataset"]]
print(dataset_curr)

ens_col_use <- seu_dataset@misc$ensembl_column_use
print(ens_col_use)

# the downtream code will not work if the seurat features aren't human ens IDs 
stopifnot(ens_col_use == "ENSG_ID")

# I will keep all thresholds although some contain >1000s of genes, to compute
# them anyways but maybe not visualize them since they are too different from
# the prupose of a short signature at first glance.

#-------------------------------------------------------------------------------

# only using signature genes
# but because the gene lists are longer and this is a human dataset
# I need the original human ensembl df as well as the mouse ensembl
# df because the human one only has the ensmus ID
ensembl_hum <- readRDS(snakemake@input[["ensembl_hum"]])
ensembl_mus <- readRDS(snakemake@input[["ensembl_mus"]])
#ensembl_hum <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/ensembl/ensembl_hum")
#ensembl_mus <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/ensembl/ensembl_mus")

#-------------------------------------------------------------------------------
# get the selected resolution

# data frame containing the correct resolutions for other datasets as visually determined
resolution_df_path <- base::paste0(
  snakemake@input[["metadata_path"]],
  "/scRNAseq/03_sce_analysis/reclustering_bm/reclustering_other_resolution.txt")

print(resolution_df_path)
resolution_df <- utils::read.csv(
  file = resolution_df_path, 
  header = TRUE, 
  sep = ";", 
  check.names=FALSE, 
  stringsAsFactors=FALSE, 
  as.is=TRUE, 
  colClasses = "character")

res_sign <- resolution_df %>%
  dplyr::filter(dataset == dataset_curr) %>%
  dplyr::filter(conservation_level == "conserved_signature") %>%
  dplyr::pull(resolution)

print(res_sign)

#-------------------------------------------------------------------------------

nr_cores <- snakemake@params[["nr_cores"]]
# nr_cores <- 2

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## for each list item of signature genes, get the corresponding ENSG IDs instead

signature_genes_by_t_ensids <- lapply(signature_genes_by_t, function(vec){
  
  print(paste("started from", length(vec)))

  musids <- ensembl_mus %>%
    dplyr::filter(external_gene_name %in% vec) %>%
    dplyr::pull(ensembl_gene_id)

  ensids <- ensembl_hum %>%
    dplyr::filter(mmusculus_homolog_ensembl_gene %in% musids) %>%
    dplyr::pull(ensembl_gene_id)
  print(paste("ended with", length(ensids)))
  
  # some genes are lost due to lack of orthologs

  return(ensids)
})

signature_genes_by_t_ensids[[1]]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# SUBSET
# get a list, each item is the same seu object subsetted by the different genes

subset_seu_list_by_t <- lapply(
  signature_genes_by_t_ensids, 
  function(ensids, seu){

    print(length(ensids))
    
    seu_sign <- BiocGenerics::subset(
      seu, 
      features = ensids,
      slot = "count")
    
    # add info as usual
    seu_sign@misc$all_features_reclustering <- ensids
    seu_sign@misc$used_genes <- "conserved_signature_treshold_testing_by_t"
    seu_sign@misc$nr_genes_used <- length(ensids)
    
    print(dim(seu_sign))
    # some more genes are lost that are not part of the object
    return(seu_sign)
  },
  seu_dataset
)

subset_seu_list_by_t[[1]]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# RECLUSTER
# based on 02_scripts_other/03_reclustering_other.R 

# go through *standard* Seurat pipeline 
# starts from raw counts 

print("nr_cores")
print(nr_cores)

# run the standard pipeline for each object, subsetted by signature genes
# obtained from different p value thresholds on nDEG analysis    
              
seu_sign_list_reclustered <- mclapply(
  X = subset_seu_list_by_t, 

  FUN = function(seu_obj){

    print(seu_obj)
    dim(seu_obj)
    seu_return <- standard_seu_pipeline( # own function
      resolution = res_sign,
      seu = seu_obj,
      data_use = seu_obj@misc$data_use,
      features = rownames(seu_obj)
    )
    
    return(seu_return)
    },
    mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE
)


# for testing
#print("testing")
#seu_sign_list_reclustered <- lapply(
#  X = subset_seu_list_by_t, 

#  FUN = function(seu_obj){

#    print(seu_obj)
#    dim(seu_obj)

#    seu_return <- standard_seu_pipeline(
#      resolution = res_sign,
#      seu = seu_obj,
#      data_use = seu_obj@misc$data_use,
#      features = rownames(seu_obj)
#    )
    
#    return(seu_return)
#  }
#)
#print("done testing")

names(seu_sign_list_reclustered)
head(seu_sign_list_reclustered[[1]]@meta.data)
seu_sign_list_reclustered[[1]]@misc

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Save

print(snakemake@output[["seu_recl_output"]])
saveRDS(seu_sign_list_reclustered, file = snakemake@output[["seu_recl_output"]])

utils::sessionInfo()