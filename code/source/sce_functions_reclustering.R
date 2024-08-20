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
    npcs = 30,
    verbose = FALSE, 
    features = features)
  
  if(ncol(seu@reductions$pca) == 30){
    nr_pca <- 30
  }else if(ncol(seu@reductions$pca) < 30){
    nr_pca <- ncol(seu@reductions$pca)
    print(base::paste("changed nr of pcr to maximum", nr_pca, "PCs"))
  }
  
  seu <- Seurat::FindNeighbors(
    seu, 
    dims = 1:nr_pca, 
    verbose = FALSE)
  
  seu <- Seurat::FindClusters(
    seu, 
    resolution = resolution, 
    verbose = FALSE) 
  
  if(calc_umap == TRUE){
    seu <- Seurat::RunUMAP(
      seu,
      dims = 1:nr_pca,
      verbose = FALSE)
  }
  
  seu@misc$resolution <- resolution
  
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
  mat <- base::table(seu$cell_type, seu$seurat_clusters)
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


calculate_scores_long <- function(
    seu = NULL, 
    cluster_vector1 = NULL,
    cluster_vector2 = NULL,
    mat_pca = NULL,
    for_permutation = FALSE){
  
  # can be seurat object, or two clustering vectors + pca matrix
  
  if(!is.null(seu)){
    seu <- seu 
    # seurat objects with clusters in "seurat_clusters" slot and
    # cell type or identity to compare in "cell_type" slot
  
    cluster_vector1 <- seu$cell_type
    cluster_vector2 <- seu$seurat_clusters
    
    mat_pca <- seu@reductions$pca@cell.embeddings[,1:10]
  }else{
    cluster_vector1 <- cluster_vector1 # old clustering
    cluster_vector2 <- cluster_vector2 # new clustering
    # purity will be calculated from cluster_vector2
  
    mat_pca <- mat_pca
  }
  
  # make a comparison matrix, then split into cells/celltype and cells/cluster
  mat <- base::table(cluster_vector1, cluster_vector2)
  mat_per_celltype <- mat/Matrix::rowSums(mat)
  t_mat <- t(mat)
  mat_per_cluster <- t(t_mat/Matrix::rowSums(t_mat)) # keep format
  
  nr_celltypes <- nrow(mat)
  nr_clusters <- ncol(mat)
  
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
  stopifnot(score_1 <= 1)
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # score 2: mean proportion of cells of a cell type per cluster
  
  mat_per_cluster[mat_per_cluster == 0] <- NA
  per_cluster_mean <- base::mean(mat_per_cluster, na.rm = TRUE)
  
  score_2 <- per_cluster_mean

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # score 3: mean purity of new clusters
  
  res_recl <- bluster::neighborPurity(mat_pca, 
                                      clusters = cluster_vector2, 
                                      k = 50)
  score_3 <- base::mean(res_recl$purity)
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # scores 4 - 6: established metrics for comparing two clusterings
  
  # Adjusted Rand index
  score_4 <- mclust::adjustedRandIndex(cluster_vector1, cluster_vector2)
  
  # Fowlkes-Mallows Index
  score_5 <- dendextend::FM_index_R(cluster_vector1, cluster_vector2)[1]
  
  # Variation of Information
  score_6 <- mcclust::vi.dist(unfactor(cluster_vector1),
                              unfactor(cluster_vector2))

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  
  # depending on what this function is used for, output differs slightly
  if(!for_permutation){
    
    return_df <- base::data.frame(
      "type" = c("proportion_of_zeros", 
                 "mean_prop_cells/cluster",
                 "mean_cluster_purity",
                 "adjusted_rand_index",
                 "fowles_mallow_index",
                 "variation_information",
                 "nr_clusters",
                 "nr_celltypes"),
      "value" = c(score_1, 
                  score_2,
                  score_3,
                  score_4, 
                  score_5,
                  score_6,
                  nr_clusters, 
                  nr_celltypes))
    
  }else if(for_permutation){
    
    return_df <- base::data.frame(
      "proportion_of_zeros" = score_1,
      "mean_prop_cells/cluster" = score_2,
      "mean_cluster_purity" = score_3,
      "adjusted_rand_index" = score_4,
      "fowles_mallow_index" = score_5,
      "variation_information" = score_6,
      "nr_clusters" = nr_clusters,
      "nr_celltypes" = nr_celltypes) 
  }
  
  return(return_df)
}  

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
# function that combines two previous functions to go through i
# iterations of seurat re-clustering using i sets of random genes and
# calculating the i sets of resulting scores
# the function will only output the scores for direct comparison 

random_reclustering_scores <- function(iteration, 
                                       seu, 
                                       iteration_df, 
                                       resolution,
                                       assay_use){
  
  iteration <- iteration # list of iterations like: as.list(c(1:i))
  seu <- seu # seurat object to be reclustered 
  iteration_df <- iteration_df # dataframe with ncol=i and random sets of gene positions in each column
  resolution <- resolution # resolution to parse to FindClusters
  assay_use <- assay_use # the assay to be used for reclustering ("RNA" in my case)
  
  # subset object to previously generated random list of genes from iteration_df
  iteration_vector <- iteration_df[,iteration]
  #print(iteration_vector[1:10])
  
  seu_sub <- SeuratObject::subset(
    seu,
    features = SeuratObject::Features(seu)[iteration_vector], 
    slot = "count")
  
  # re-cluster seurat object from raw counts, using basic seurat pipeline
  seu_rec <- standard_seu_pipeline(
    seu = seu_sub, 
    features = SeuratObject::Features(seu_sub), 
    resolution = resolution,
    assay_use = assay_use)
  
  # calculate re-clustering scores as defined above
  res_df <- calculate_scores_long(seu_rec)
  res_df$iteration[1] <-iteration
  
  return(res_df)
}



