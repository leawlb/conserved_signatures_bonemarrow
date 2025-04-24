# get reclustering stores for reclustered datasets from other species

# determine random number generator for sample
library(parallel)
RNGkind("L'Ecuyer-CMRG") # using this for reclustering

fraction_curr <- snakemake@wildcards[["fraction"]]

#-------------------------------------------------------------------------------

library(mclust, quietly = TRUE)
library(bluster, quietly = TRUE)
library(S4Vectors, quietly = TRUE)

source(file = snakemake@params[["functions_reclustering"]])

#-------------------------------------------------------------------------------
# load seurat object lists with reclustered labels

sce_list_reso <- base::readRDS(snakemake@input[["sce_list_reso"]])

#-------------------------------------------------------------------------------

# add info to each sce object

for(res in names(sce_list_reso)){
  print(res)
  for(cons_level in names(sce_list_reso[[res]])){
    
    #print(cons_level)
    #print(sce_list_reso[[res]][[cons_level]])
    
    if(cons_level == "sce_signt"){
      cl <- "conserved_signature"
    }else if(cons_level == "sce_consm"){
      cl <- "conserved_markers"
    }else if(cons_level == "sce_mmusm"){
      cl <- "mmusall_markers"
    }
    
    sce_list_reso[[res]][[cons_level]]$conservation_level <- cl
    sce_list_reso[[res]][[cons_level]]$resolution <- res
  
  }
}
print(sce_list_reso[[1]][[1]]$conservation_level[1])
print(sce_list_reso[[1]][[1]]$resolution[1])
print(sce_list_reso[[1]][[1]]$nr_genes[1])

#-------------------------------------------------------------------------------
# for each item of the nested list, calculate the reclustering scores 
score_df_list_all <- lapply(sce_list_reso, function(sce_list){
  
  score_df_list <- lapply(sce_list, function(sce){
    
    cluster_vector1 <- sce$celltypes # cluster_vector1 = cell types
    cluster_vector2 <- sce$reclustered # cluster_vector2 = new clusters
    mat_pca <- SingleCellExperiment::reducedDims(sce)$PCA[,1:10]
    
    # use own function for calculation
    score_df <- calculate_scores( 
      cluster_vector1 = cluster_vector1,
      cluster_vector2 = cluster_vector2,
      mat_pca = mat_pca)
    
    # add specific info
    score_df$conservation_level <- base::rep(sce$conservation_level[1],
                                             nrow(score_df))
    score_df$resolution <- base::rep(sce$resolution[1],
                                     nrow(score_df))
    score_df$nr_genes_used <- base::rep(sce$nr_genes[1],
                                     nrow(score_df))
    
    score_df$nr_celltypes <- base::rep(
      length(base::unique(cluster_vector1)), 
      nrow(score_df))
    score_df$nr_clusters <- base::rep(
      length(base::unique(cluster_vector2)), 
      nrow(score_df))
    score_df$fraction <- base::rep(
      fraction_curr, 
      nrow(score_df))

    
    return(score_df)
  })
  names(score_df_list) <- names(sce_list)
  return(score_df_list)
})
names(score_df_list_all) <- names(sce_list_reso)

#-------------------------------------------------------------------------------

base::saveRDS(score_df_list_all, snakemake@output[["score_df_list"]])

utils::sessionInfo()
