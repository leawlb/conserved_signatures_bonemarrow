#-------------------------------------------------------------------------------
# re-cluster own datasets with each gene set 
# - conserved signature,
# - conserved markers
# - all BL6 genes,
# - nDGEs

# subclustering genes are removed in cell-type specific manner

#-------------------------------------------------------------------------------

# determine random number generator for sample
# Mersenne-Twister" is default
RNGkind("Mersenne-Twister") 

set.seed(37)

library(scater, quietly = TRUE)
library(scran, quietly = TRUE)
library(bluster, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(parallel, quietly = TRUE)

source(snakemake@params[["functions_reclustering"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# load

# sce object
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

# genesets
geneset_list <- base::readRDS(snakemake@input[["geneset_list"]])

#-------------------------------------------------------------------------------
# params for reclustering
k_graph_list <- snakemake@params[["k_graph_list"]]
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]
resolution_louvain <- resolution_louvain_list[[fraction_curr]]

nr_cores <- snakemake@params[["nr_cores"]]

print(k_graph)
print(resolution_louvain)
print(nr_cores)

#-------------------------------------------------------------------------------
# cell types that were not used to extract signature are excluded because
# they cannot be separated after excluding them
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE

#-------------------------------------------------------------------------------
# get the required gene sets
# remove subclustering genes from gene sets in a cell-type specific manner
# this allows for genes that were used to sub-cluster two cell types to
# still be used for re-clustering if they are part of the conserved signature
# of another cell type

# get all signature genes
conserved_signature_list <- lapply(geneset_list, function(geneset){
  conserved_signature <- geneset$conserved_signature
  conserved_signature <- conserved_signature[
    which(!conserved_signature %in% geneset$subclustering_genes)]
  return(conserved_signature)
})
signt <- base::unique(unlist(conserved_signature_list))
print(length(signt))

# get all conserved markers
conserved_markers_list <- lapply(geneset_list, function(geneset){
  conserved_markers <- geneset$conserved_markers
  conserved_markers <- conserved_markers[
    which(!conserved_markers %in% geneset$subclustering_genes)]
  return(conserved_markers)
})
consm <- base::unique(unlist(conserved_markers_list))
print(length(consm))

# get all ndges
ndge_list <- lapply(geneset_list, function(geneset){
  ndges <- geneset$ndges
  ndges <- ndges[
    which(!ndges %in% geneset$subclustering_genes)]
  return(ndges)
})
ndges <- base::unique(unlist(ndge_list))
print(length(ndges))

# get all BL6 markers
mmus_marker_list <- lapply(geneset_list, function(geneset){
  mmus_markers <- geneset$conserved_df$gene[
    which(!is.na(geneset$conserved_df$mmus))]
  mmus_markers <- mmus_markers[
    which(!mmus_markers %in% geneset$subclustering_genes)]
  return(mmus_markers)
})
mmusm <- base::unique(unlist(mmus_marker_list))
print(length(mmusm))

# extract subclustering genes
subclustering_gene_list <- lapply(geneset_list, function(geneset){
  subclustering_genes <- geneset$genes_subclustering
  return(subclustering_genes)
})
subclustering_genes <- base::unique(unlist(subclustering_gene_list))

#-------------------------------------------------------------------------------
# subset to required cell types
sce <- sce[,which(!sce$celltypes %in% cts_exclude)]

# remove cell type levels that are not required anymore
sce$celltypes <- factor(
  sce$celltypes,
  levels = levels(sce$celltypes)[
    which(levels(sce$celltypes) %in% base::unique(sce$celltypes))])

# if fraction is hsc, downsample to accelerate re-clustering
if(fraction_curr == "hsc"){
  # set a different seed so the same numbers are not generated accidentally 
  set.seed(147)
  
  downample_pos <- base::sample(c(1:ncol(sce)), 25000, replace = FALSE)
  sce <- sce[,downample_pos]
  print(dim(sce))
  print(head(downample_pos))
  
  # set seed back to seed
  set.seed(37)
}

#-------------------------------------------------------------------------------

# subset to all genes that are not subclustering genes
# TODO: if ok, remove the cell type-specific removal
sce <- sce[which(!rownames(sce) %in% subclustering_genes)]
print(dim(sce))

# subset to gene sets
sce_signt <- sce[rownames(sce) %in% signt,]
sce_consm <- sce[rownames(sce) %in% consm,]
sce_ndges <- sce[rownames(sce) %in% ndges,]
sce_mmusm <- sce[rownames(sce) %in% mmusm,]

print(dim(sce_signt))
print(dim(sce_consm))
print(dim(sce_ndges))
print(dim(sce_mmusm))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# RE-CLUSTER

#-------------------------------------------------------------------------------
# re-cluster and add clusters to SCE object
# exactly how the data was orignally clustered

# cluster each subsetted sce object
sce_list <- list("sce_signt" = sce_signt,
                 "sce_consm" = sce_consm, 
                 "sce_ndges" = sce_ndges,
                 "sce_mmusm" = sce_mmusm)

clustered_list <- mclapply(
  sce_list,
  clustering_orig,
  k_graph = k_graph, 
  resolution_louvain = resolution_louvain,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE)

# for testing without mclapply
# clustered_list <- lapply(
#   sce_list,
#   clustering_orig,
#   k_graph = k_graph, 
#   resolution_louvain = resolution_louvain)

#-------------------------------------------------------------------------------
# transfer clustering info to original sce

sce$cluster_consm <- clustered_list$sce_consm$reclustered
sce$cluster_ndges <- clustered_list$sce_ndges$reclustered
sce$cluster_signt <- clustered_list$sce_signt$reclustered
sce$cluster_mmusm <- clustered_list$sce_mmusm$reclustered

print(head(colData(sce)))

#-------------------------------------------------------------------------------
base::saveRDS(sce, file = snakemake@output[["sce_output"]])

utils::sessionInfo()
