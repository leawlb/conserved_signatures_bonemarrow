#-------------------------------------------------------------------------------
# re-cluster own datasets with each gene set 
# - conserved signature,
# - conserved markers,
# - all BL6 genes,
# - nDGEs

# subclustering genes are removed IN TOTAL

#-------------------------------------------------------------------------------

# determine random number generator for sample
library(parallel)
RNGkind("L'Ecuyer-CMRG") # using this for parallel is necessary

set.seed(37)

library(scater, quietly = TRUE)
library(scran, quietly = TRUE)
library(bluster, quietly = TRUE)
library(tidyverse, quietly = TRUE)

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
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]

nr_cores <- snakemake@params[["nr_cores"]]

print(k_graph)
print(nr_cores)

#-------------------------------------------------------------------------------
# cell types that were not used to extract signature are excluded
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE

#-------------------------------------------------------------------------------
# get the required gene sets
# sub-clustering genes will be removed IN TOTAL downstream

# get all signature genes
conserved_signature_list <- lapply(geneset_list, function(geneset){
  conserved_signature <- geneset$conserved_signature
  return(conserved_signature)
})
signt <- base::unique(unlist(conserved_signature_list))
print(length(signt))

# get all conserved markers
conserved_markers_list <- lapply(geneset_list, function(geneset){
  conserved_markers <- geneset$conserved_markers
  return(conserved_markers)
})
consm <- base::unique(unlist(conserved_markers_list))
print(length(consm))

# get all ndges
ndge_list <- lapply(geneset_list, function(geneset){
  ndges <- geneset$ndges
  return(ndges)
})
ndges <- base::unique(unlist(ndge_list))
print(length(ndges))

# get all BL6 markers
mmus_marker_list <- lapply(geneset_list, function(geneset){
  mmus_markers <- geneset$conserved_df$gene[
    which(!is.na(geneset$conserved_df$mmus))]
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
  set.seed(14744)
  
  downsample_pos <- base::sample(c(1:ncol(sce)), 25000, replace = FALSE)
  sce <- sce[,downsample_pos]
  print(dim(sce))
  print(head(downsample_pos))
  
  # set seed back to seed
  set.seed(37)
}

#-------------------------------------------------------------------------------

# subset to all genes that are not subclustering genes
# SUBCLUSTERING GENES REMOVED IN TOTAL
print(subclustering_genes)
print(length(subclustering_genes))
sce <- sce[which(!rownames(sce) %in% subclustering_genes)]
print(dim(sce))

# remove subclustering genes from the object
# this will also ensure that all downstream re-clusterings based on this 
# exported objects will have sub-clustering genes removed
sce_signt <- sce[rownames(sce) %in% signt,]
sce_consm <- sce[rownames(sce) %in% consm,]
sce_mmusm <- sce[rownames(sce) %in% mmusm,]

print(dim(sce_signt))
print(dim(sce_consm))
print(dim(sce_mmusm))

stopifnot(!rownames(sce_signt) %in% subclustering_genes)
stopifnot(!rownames(sce_consm) %in% subclustering_genes)
stopifnot(!rownames(sce_mmusm) %in% subclustering_genes)

sce_signt$nr_genes <- length(signt)
sce_consm$nr_genes <- length(consm)
sce_mmusm$nr_genes <- length(mmusm)
print(sce_signt$nr_genes[1])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# RE-CLUSTER

# so far exactly as 03_reclustering_own.R, except removal of nDEGs
# now, introduce resolutions_list

resolution_list <- as.list(seq(0.1, 1.3, by = 0.05))

print("nr_cores")
print(nr_cores)

set.seed(37)

# cluster each subsetted sce object
# leave out nDGEs
sce_list <- list("sce_signt" = sce_signt,
                 "sce_consm" = sce_consm, 
                 "sce_mmusm" = sce_mmusm)

reso_reclustered_list <- mclapply(resolution_list, function(resolution){
  
  res <- resolution
  print(res)
  print(k_graph)
  print(sce_list)
  
  reso_recl <- lapply(
    sce_list, # four sce objects
    clustering_orig, # own function
    k_graph = k_graph, 
    resolution_louvain = res)
  
  names(reso_recl) <- names(sce_list)
  
  return(reso_recl)
  
  },
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE
  )

names(reso_reclustered_list) <- as.character(unlist(resolution_list))

base::saveRDS(reso_reclustered_list, 
              snakemake@output[["reso_reclustered_list"]])

utils::sessionInfo()
