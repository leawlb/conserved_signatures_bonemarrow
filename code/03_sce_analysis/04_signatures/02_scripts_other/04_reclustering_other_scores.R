
# get reclustering stores for reclustered datasets from other species

# determine random number generator for sample
library(parallel)
RNGkind("L'Ecuyer-CMRG") # using this for reclustering

#-------------------------------------------------------------------------------

library(mclust, quietly = TRUE)
library(bluster, quietly = TRUE)
library(S4Vectors, quietly = TRUE)

source(file = snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
# load seurat object lists with reclustered labels

seu_list_all <- base::readRDS(snakemake@input[["seu_list"]])

#-------------------------------------------------------------------------------
# for each item of the nested list, calculate the reclustering scores 
score_df_list_all <- lapply(seu_list_all, function(seu_list){
  
  score_df_list <- lapply(seu_list, function(seu){
    
    print(seu@misc$used_genes)
    print(seu@misc$resolution)
    score_df <- calculate_scores(seu) # own function
    
    # add specific info
    score_df$conservation_level <- base::rep(seu@misc$used_genes,
                                             nrow(score_df))
    score_df$resolution <- base::rep(seu@misc$resolution,
                                             nrow(score_df))   
    score_df$nr_genes_used <- base::rep(seu@misc$nr_genes_used,
                                        nrow(score_df))
    
    return(score_df)
  })
  names(score_df_list) <- names(seu_list)
  return(score_df_list)
})
names(score_df_list_all) <- names(seu_list_all)

#-------------------------------------------------------------------------------

base::saveRDS(score_df_list_all, snakemake@output[["score_df_list"]])

utils::sessionInfo()
