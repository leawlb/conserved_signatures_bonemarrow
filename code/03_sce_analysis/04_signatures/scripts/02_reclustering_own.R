#-------------------------------------------------------------------------------
# re-cluster our own two datasets with each gene set (cons markers, nDGEs, EMFs)

#-------------------------------------------------------------------------------

library(scater, quietly = TRUE)
library(scran, quietly = TRUE)
library(bluster, quietly = TRUE)
library(tidyverse, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# load

signature_list <- base::readRDS(snakemake@input[["signature_list"]])
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

k_graph_list <- snakemake@params[["k_graph_list"]]
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]
resolution_louvain <- resolution_louvain_list[[fraction_curr]]
print(k_graph)
print(resolution_louvain)

# cell types that were not used to extract signature are excluded because
# they cannot be separated after excluding them
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)

#-------------------------------------------------------------------------------
# get the required gene sets
# remove subclustering genes from gene sets in a cell-type specific manner
# this allows for genes that were used to sub-cluster two cell types to
# still be used for re-clustering if they are part of the conserved signature
# of another cell type

# get all conserved markers
conserved_markers_list <- lapply(signature_list, function(sign){
  conserved_markers <- sign$conserved_markers
  conserved_markers <- conserved_markers[
    !conserved_markers %in% sign$subclustering_genes]
  return(conserved_markers)
})
consm <- base::unique(unlist(conserved_markers_list))

# get all ndges
ndge_list <- lapply(signature_list, function(sign){
  ndges <- sign$ndges
  ndges <- ndges[!ndges %in% sign$subclustering_genes]
  return(ndges)
})
ndges <- base::unique(unlist(ndge_list))

# get all signature genes
conserved_signature_list <- lapply(signature_list, function(sign){
  conserved_signature <- sign$conserved_signature
  conserved_signature <- conserved_signature[
    !conserved_signature %in% sign$subclustering_genes]
  return(conserved_signature)
})
signt <- base::unique(unlist(conserved_signature_list))

#-------------------------------------------------------------------------------
# subset to required cell types
sce <- sce[,!sce$celltypes %in% cts_exclude]

# subset to gene sets
sce_consm <- sce[rownames(sce) %in% consm,]
sce_ndges <- sce[rownames(sce) %in% ndges,]
sce_signt <- sce[rownames(sce) %in% signt,]

print(dim(sce_consm))
print(dim(sce_ndges))
print(dim(sce_signt))

#-------------------------------------------------------------------------------
# re-cluster and add clusters to SCE object
# exactly how the data was orignally clustered

# make functions
clustering_orig <- function(sce, k_graph, resolution_louvain){
  
  sce <- sce
  k_graph <- k_graph
  resolution_louvain <- resolution_louvain
  
  print(nrow(sce))
  # run PCA 
  # re-calculate PCA (without batch correction)
  # use logcounts assay, contains normalized logcounts from multibatchnorm
  sce <- scater::runPCA(sce, 
                        ncomponents = 25, 
                        exprs_values = "logcounts") 
  
  print(nrow(sce))
  # get a graph of nearest neighbors
  graph <- scran::buildSNNGraph(sce, 
                                k = k_graph, 
                                use.dimred = "PCA", 
                                type = "rank")
  
  # community detection
  clust <- igraph::cluster_louvain(graph, resolution = resolution_louvain)
  
  # add to sce object 
  sce$reclustered <- clust$membership
  
  sce$reclustered <- base::factor(
    sce$reclustered, 
    levels = base::sort(base::unique(sce$reclustered)))
  print(base::table(sce$reclustered))
  
  sce$k_graph <- base::rep(k_graph, ncol(sce))
  sce$resolution_louvain <- base::rep(resolution_louvain, ncol(sce))
  
  return(sce)
}

#-------------------------------------------------------------------------------

# cluster each subset sce object
sce_list <- list("sce_consm" = sce_consm, 
                 "sce_ndges" = sce_ndges,
                 "sce_signt" = sce_signt)
clustered_list <- lapply(sce_list, function(sce, k_graph, resolution_louvain){
  
  sce <- sce
  k_graph <- k_graph
  resolution_louvain <- resolution_louvain
  
  print(nrow(sce))
  sce_clust <- clustering_orig(sce, 
                               k_graph = k_graph,
                               resolution_louvain = resolution_louvain)
  return(sce_clust)
  
}, k_graph = k_graph, resolution_louvain = resolution_louvain)

#-------------------------------------------------------------------------------
# transfer clustering info to orignal sce

sce$cluster_consm <- clustered_list$sce_consm$reclustered
sce$cluster_ndges <- clustered_list$sce_ndges$reclustered
sce$cluster_signt <- clustered_list$sce_signt$reclustered

print(head(colData(sce)))
base::saveRDS(sce, file = snakemake@output[["sce_output"]])

utils::sessionInfo()
