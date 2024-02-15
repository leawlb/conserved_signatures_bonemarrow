#-------------------------------------------------------------------------------

library(scater)
library(scran)
library(bluster)
library(tidyverse)
set.seed(37)

#-------------------------------------------------------------------------------

core_cons_list <- readRDS(snakemake@input[["core_cons_list"]])
sce <- readRDS(file = snakemake@input[["sce_input"]])

k_louvain <- snakemake@params[["k_louvain"]]

fraction_curr <- snakemake@wildcards[["fraction"]]

core_cons_df_list <- lapply(as.list(names(core_cons_list)), function(c){
  
  core_cons_list[[c]]$marker_cons_ct$cell_type <- c
  if(fraction_curr == "hsc"){
    core_cons_list[[c]]$marker_cons_ct$assignment <- "receiver"
  }else if(fraction_curr == "str"){
    core_cons_list[[c]]$marker_cons_ct$assignment <- "emitter"
  }
  
  return(core_cons_list[[c]]$marker_cons_ct)
  
})
ndge_list <- lapply(as.list(names(core_cons_list)), function(c){
  
  ndges <- core_cons_list[[c]]$ndges
  return(ndges)
  
})
ndges <- unique(unlist(ndge_list))

core_cons_df <- bind_rows(core_cons_df_list)

cc_list <- unique(core_cons_df$gene[which(core_cons_df$cons_emf == "TRUE")])
nd_list <- ndges
cm_list <- unique(core_cons_df$gene[which(core_cons_df$cons_marker == "TRUE")])

#-------------------------------------------------------------------------------

# subset to only relevant genes
sce_cc <- sce[rownames(sce) %in% cc_list,]
sce_nd <- sce[rownames(sce) %in% nd_list,]
sce_cm <- sce[rownames(sce) %in% cm_list,]

#-------------------------------------------------------------------------------

# re-calculate PCA (Without batch correction)
sce_cc <- runPCA(sce_cc, ncomponents=25) 
sce_nd <- runPCA(sce_nd, ncomponents=25) 
sce_cm <- runPCA(sce_cm, ncomponents=25) 

# calculate clusters
louvain_cc <- clusterCells(sce_cc, use.dimred="PCA", assay.type = NULL,
                           BLUSPARAM=NNGraphParam(cluster.fun="louvain", 
                                                  k = k_louvain))

louvain_nd <- clusterCells(sce_nd, use.dimred="PCA", assay.type = NULL,
                           BLUSPARAM=NNGraphParam(cluster.fun="louvain", 
                                                  k = k_louvain))

louvain_cm <- clusterCells(sce_cm, use.dimred="PCA", assay.type = NULL,
                           BLUSPARAM=NNGraphParam(cluster.fun="louvain", 
                                                  k = k_louvain))

sce$cluster_core <- louvain_cc
sce$cluster_ndge <- louvain_nd
sce$cluster_cmrk <- louvain_cm

sce$cluster_core <- factor(sce$cluster_core,
                           levels = sort(unique(sce$cluster_core))) 
sce$cluster_ndge <- factor(sce$cluster_ndge,
                           levels = sort(unique(sce$cluster_ndge))) 
sce$cluster_cmrk <- factor(sce$cluster_cmrk,
                           levels = sort(unique(sce$cluster_cmrk))) 

sce$nr_cc <- nrow(sce_cc)
sce$nr_nd <- nrow(sce_nd)
sce$nr_cm <- nrow(sce_cm)

saveRDS(sce, file = snakemake@output[["sce_output"]])