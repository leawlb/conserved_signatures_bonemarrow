#-------------------------------------------------------------------------------
# function for a standard seurat pipeline clustering pipeline from counts
# make it so that mapply can be used for list of different resolution(s)
# used for reclustering, and for reclustering permutation
standard_seu_pipeline <- function(resolution, 
                                  features, 
                                  seu, 
                                  assay_use, # e.g. "RNA"
                                  use_raw_counts = TRUE,
                                  calc_umap = FALSE){
  
  seu <- seu # seurat object
  features <- features # features for re-clustering (all Features(seu))
  resolution <- resolution # resolution for FindClusters
  assay_use <- assay_use # assay to use as basis for clustering
  calc_umap <- calc_umap # whether umap coordinates should be calculated
  
  # NormalizeData normalizes "count data"
  seu <- Seurat::NormalizeData(
    seu,
    assay = assay_use, 
    verbose = FALSE)
  
  seu <- Seurat::ScaleData(
    seu, 
    verbose = FALSE)
  
  seu <- Seurat::RunPCA(
    seu,
    verbose = FALSE, 
    features = features)
  
  seu <- Seurat::FindNeighbors(
    seu, 
    dims = 1:30, 
    verbose = FALSE)
  
  seu <- Seurat::FindClusters(
    seu, 
    resolution = resolution, 
    verbose = FALSE) 
  
  if(calc_umap == TRUE){
    seu <- Seurat::RunUMAP(
      seu,
      dims = 1:30,
      verbose = FALSE)
  }
  
  return(seu)
}

#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# function for calculating reclustering "scores"

calculate_scores <- function(seu){
  
  seu <- seu 
  # seurat objects with clusters in "seurat_clusters" slot and
  # cell type or identity to compare in "cell_type" slot
  
  # make a comparison matrix, then split into cells/celltype and cells/cluster
  mat <- table(seu$cell_type, seu$seurat_clusters)
  mat_per_celltype <- mat/rowSums(mat)
  t_mat <- t(mat)
  mat_per_cluster <- t(t_mat/rowSums(t_mat)) # keep format/orientation
  
  #-----------------------------------------------------------------------------
  # score 1: proportions of 0
  
  # extract which nr is higher; nr of clusters or number of cell types  
  if(nrow(mat) >= ncol(mat)){
    max_nr <- nrow(mat)
  }else if(nrow(mat) < ncol(mat)){
    max_nr <- ncol(mat)
  }
  
  # get the proportion of fields that are 0 from the max possible nr of fields
  # that can be 0 in a theoretical perfect re-clustering. 
  # That "max possible nr" is equivalent to all fields - max_nr
  
  score_1 <- length(which(mat == 0))/(ncol(mat)*nrow(mat)-max_nr)
  
  #-----------------------------------------------------------------------------
  # score 2: mean proportion of cells/cell type and cells/cluster
  
  # for each field
  mat_mean_prop <- (mat_per_celltype+mat_per_cluster)/2
  
  # average across all fields
  non_zeros <- mat_mean_prop[which(mat_mean_prop > 0)]
  score_2 <- sum(non_zeros)/length(non_zeros)
  
  return(
    data.frame(
      "score_1" = score_1,
      "score_2" = score_2))
}  

#-----------------------------------------------------------------------------




#-------------------------------------------------------------------------------
# function that combines the two previous functions to go through i
# iterations of seurat re-clustering using i sets of random genes and
# calculating the i sets of resulting scores
# the function will only output the scores for direct comparison 

random_reclustering_scores <- function(iteration, 
                                       seu, 
                                       iteration_df, 
                                       resolution,
                                       assay_use){
  
  iteration <- iteration # list of iterations like: as.list(c(1:i))
  #print(iteration)
  seu <- seu # seurat object to be reclustered 
  iteration_df <- iteration_df # dataframe with ncol=i and random sets of gene positions in each column
  resolution <- resolution # resolution to parse to FindClusters
  assay_use <- assay_use # the assay to be used for reclustering ("RNA" in my case)
  
  # subset object to previously generated random list of genes from iteration_df
  iteration_vector <- iteration_df[,iteration]
  #print(iteration_vector[1:10])
  
  seu_sub <- subset(seu, features = Features(seu)[iteration_vector], 
                    slot = "count")
  
  # re-cluster seurat object from raw counts, using basic seurat pipeline
  seu_rec <- standard_seu_pipeline(seu = seu_sub, 
                                   features = Features(seu_sub), 
                                   resolution = resolution,
                                   assay_use = assay_use)
  
  # calculate re-clustering scores as defined above
  res_df <- calculate_scores(seu_rec)
  res_df$iteration[1] <-iteration
  
  return(res_df)
}

