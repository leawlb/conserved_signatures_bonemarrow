#-------------------------------------------------------------------------------

# from https://www.bioconductor.org/packages/devel/bioc/vignettes/scmap/inst/doc/scmap.html
# uses logcounts (not batch corrected) values as far as I'm aware

library(DropletUtils, quietly = TRUE) 
library(scmap, quietly = TRUE)
set.seed(37)
#-------------------------------------------------------------------------------

sce_12 <- readRDS(file = snakemake@input[["sce_12"]])
sce_07 <- readRDS(file = snakemake@input[["sce_07"]])
sce_mmus <- readRDS(file = snakemake@input[["mmus"]])

print(sce_12)
print(sce_07)
print(sce_mmus)

#-------------------------------------------------------------------------------
# use sce_07 as it is not batch corrected and therefore contains all genes
# scmap doesn't use batch corrected values to my knowledge

sce_07$cluster_hierarchical <- sce_12$cluster_hierarchical[match(colnames(sce_07), colnames(sce_12))]
sce_07$cluster_seurat <- sce_12$cluster_seurat[match(colnames(sce_07), colnames(sce_12))]
sce_07$cluster_louvain <- sce_12$cluster_louvain[match(colnames(sce_07), colnames(sce_12))]

#-------------------------------------------------------------------------------
## SCMAP Cluster

# use mmus as reference for cluster label transfer
sce_ref <- sce_mmus
rowData(sce_ref)$feature_symbol <- rownames(sce_ref)
sce_ref <- selectFeatures(sce_ref, suppress_plot = TRUE)

get_scmapcluster <- function(sce){
  
  print("hierarchical")
  sce_ref_temp1 <- sce_ref # reference data
  sce_ref_temp1 <- indexCluster(sce_ref_temp1, 
                                cluster_col = "cluster_hierarchical")
  
  sce_test_temp1 <- sce # test data
  rowData(sce_test_temp1)$feature_symbol <- rownames(sce_test_temp1)
  results <- scmapCluster(projection = sce_test_temp1,
                          index_list = list(
                            metadata(sce_ref_temp1)$scmap_cluster_index), 
                          threshold = 0.7)
  
  sce$mmus_hierarchical_pred_cluster <- results$combined_labs
  sce$mmus_hierarchical_pred_cluster_sim <- results$scmap_cluster_siml
  
  print("seurat")
  sce_ref_temp2 <- sce_ref # reference data
  sce_ref_temp2 <- indexCluster(sce_ref_temp2, cluster_col = "cluster_seurat")
  
  sce_test_temp2 <- sce # test data
  rowData(sce_test_temp2)$feature_symbol <- rownames(sce_test_temp2)
  results <- scmapCluster(projection = sce_test_temp2,
                          index_list = list(
                            metadata(sce_ref_temp2)$scmap_cluster_index), 
                          threshold = 0.7)
  
  sce$mmus_seurat_pred_cluster <- results$combined_labs
  sce$mmus_seurat_pred_cluster_sim <- results$scmap_cluster_siml
  
  print("louvain")
  sce_ref_temp3 <- sce_ref # reference data
  sce_ref_temp3 <- indexCluster(sce_ref_temp3, cluster_col = "cluster_louvain")
  
  sce_test_temp3 <- sce # test data
  rowData(sce_test_temp3)$feature_symbol <- rownames(sce_test_temp3)
  results <- scmapCluster(projection = sce_test_temp3,
                          index_list = list(
                            metadata(sce_ref_temp3)$scmap_cluster_index), 
                          threshold = 0.7)
  
  sce$mmus_louvain_pred_cluster <- results$combined_labs
  sce$mmus_louvain_pred_cluster_sim <- results$scmap_cluster_siml
  
  return(sce)
}

print("starting get_scmapcluster")
sce_14_nbc <- get_scmapcluster(sce_07)
sce_14 <- get_scmapcluster(sce_12)

#-------------------------------------------------------------------------------
## SCMAP Cell

# use mmus as reference for cell projection and cluster label transfer

# stochastic
set.seed(37)
sce_ref2 <- sce_mmus
rowData(sce_ref2)$feature_symbol <- rownames(sce_ref2)
sce_ref2 <- selectFeatures(sce_ref2, suppress_plot = TRUE)
sce_ref2 <- indexCell(sce_ref2)

# in https://www.bioconductor.org/packages/devel/bioc/vignettes/scmap/inst/doc/scmap.html
# only set.seed is used, but no random numbers are used as far as I'm aware (?)
get_scmapcell <- function(sce){
  
  sce_test_temp <- sce # test data
  rowData(sce_test_temp)$feature_symbol <- rownames(sce_test_temp)
  
  # project cells
  scmapCell_results <- scmapCell(
    projection = sce_test_temp, 
    index_list = list(
      results = metadata(sce_ref2)$scmap_cell_index
    )
  )
  
  # get all cluster labels
  print("hierarchical")
  scmapCell_clusters_hir <- scmapCell2Cluster(
    scmapCell_results = scmapCell_results, 
    cluster_list = list(
      as.character(colData(sce_ref2)$cluster_hierarchical)
    )
  )
  
  print("seurat")
  scmapCell_clusters_ser <- scmapCell2Cluster(
    scmapCell_results = scmapCell_results, 
    cluster_list = list(
      as.character(colData(sce_ref2)$cluster_seurat)
    )
  )
  
  print("louvain")
  scmapCell_clusters_lou <- scmapCell2Cluster(
    scmapCell_results = scmapCell_results, 
    cluster_list = list(
      as.character(colData(sce_ref2)$cluster_louvain)
    )
  )
  
  print(scmapCell_clusters_hir$scmap_cluster_siml)
  
  # add to sce 
  sce$mmus_hierarchical_pred_cell <- scmapCell_clusters_hir$combined_labs
  sce$mmus_hierarchical_pred_cell_sim <-scmapCell_clusters_hir$scmap_cluster_siml
  
  sce$mmus_seurat_pred_cell <- scmapCell_clusters_ser$combined_labs
  sce$mmus_seurat_pred_cell_sim <- scmapCell_clusters_ser$scmap_cluster_siml
  
  sce$mmus_louvain_pred_cell <- scmapCell_clusters_lou$combined_labs
  sce$mmus_louvain_pred_cell_sim <- scmapCell_clusters_lou$scmap_cluster_siml
  
  return(sce)
}  

print("starting get_scmapcell")
sce_14_nbc <- get_scmapcell(sce_07)
sce_14 <- get_scmapcell(sce_14)

print(sce_14_nbc)
print(sce_14)
print(colData(sce_14_nbc))
print(colData(sce_14))

saveRDS(sce_14, file = snakemake@output[["sce_14"]])
saveRDS(sce_14_nbc, file = snakemake@output[["sce_14_nbc"]])
