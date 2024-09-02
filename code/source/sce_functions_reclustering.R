
# functions used for re-clustering of or permutation tests for own datasets,
# or other datasets using different gene sets

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

######################### standard_seu_pipeline ################################

#-------------------------------------------------------------------------------

# function for a standard seurat pipeline clustering pipeline from raw counts
# or, if counts are not available, starting at scaling normalised counts

# designed so that mapply can be used for a list of different resolution(s)

standard_seu_pipeline <- function(resolution, 
                                  features, 
                                  seu, 
                                  data_use, 
                                  calc_umap = FALSE){
  
  seu <- seu # seurat object
  features <- features # features for (re-)clustering 
  resolution <- resolution # resolution for FindClusters
  data_use <- data_use # assay to use as basis for clustering, either raw_counts
                       # or log_counts
  calc_umap <- calc_umap # whether umap coordinates should be calculated
  
  print(resolution)
  print(data_use)
  print("......starting Seurat pipeline")
  
  if(data_use == "raw_counts"){
    
    # NormalizeData normalizes "count data"
    seu <- Seurat::NormalizeData(
      seu,
      assay = "RNA", 
      verbose = FALSE)
  }
  
  # data needs to be scaled either way
  seu <- Seurat::ScaleData(
    seu, 
    verbose = FALSE)
  
  # uses scaled data
  seu <- Seurat::RunPCA(
    seu,
    npcs = 30,
    verbose = FALSE, 
    features = features)
  
  # in case the number of PCs is lower than 30 due to small datasets
  if(ncol(seu@reductions$pca) == 30){
    nr_pca <- 30
  }else if(ncol(seu@reductions$pca) < 30){
    nr_pca <- ncol(seu@reductions$pca)
    print(base::paste("changed nr of pcr to maximum", nr_pca, "PCs"))
  }
  print(nr_pca)
  
  # requires PC coordinates
  seu <- Seurat::FindNeighbors(
    seu, 
    reduction = "pca",
    dims = 1:nr_pca, 
    features = features,
    verbose = FALSE)
  
  # requires neighbors
  seu <- Seurat::FindClusters(
    seu, 
    resolution = as.numeric(resolution), 
    verbose = FALSE) 
  
  if(calc_umap == TRUE){
    print("calculating umap")
    seu <- Seurat::RunUMAP(
      seu,
      dims = 1:nr_pca,
      verbose = FALSE)
  }
  
  print(".....finished Seurat pipeline")
  seu@misc$resolution <- resolution
  
  return(seu)
}



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

############################## clustering_orig #################################

#-------------------------------------------------------------------------------

# original clustering pipeline as used for our datasets
# used here for re-clustering own datasets as part of permutation and 
# gene set testing

clustering_orig <- function(sce, k_graph, resolution_louvain){
  
  sce <- sce # SCE object (own dataset)
  k_graph <- k_graph # k = number of nearest neighbors for buildSNNGraph()
  resolution_louvain <- resolution_louvain # resolution for cluster_louvain()
  
  print(nrow(sce))
  print(ncol(sce))
  
  # re-calculate PCA (without batch correction)
  # use logcounts assay, contains normalized logcounts from multibatchnorm
  # remove old PCA
  reducedDims(sce)$PCA <- NULL
  sce <- scater::runPCA(sce, 
                        ncomponents = 25, 
                        exprs_values = "logcounts") 
  
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
  print(base::table(sce$reclustered, sce$celltypes))
  
  sce$k_graph <- base::rep(k_graph, ncol(sce))
  sce$resolution_louvain <- base::rep(resolution_louvain, ncol(sce))
  
  return(sce)
}



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

############################## calculate_scores ################################

#-------------------------------------------------------------------------------

# calculate chosen scores for evaluating re-clustering

calculate_scores <- function(
    seu = NULL, 
    cluster_vector1 = NULL, 
    cluster_vector2 = NULL, 
    mat_pca = NULL, 
    for_permutation = FALSE){ 
  
  seu <- seu # seurat object to be evaluated
  cluster_vector1 <- cluster_vector1 # or, a vector of original clusters 
  cluster_vector2 <- cluster_vector2 # and of new clusters for comparison
  mat_pca <- mat_pca # and a matrix of PC coordinates
  for_permutation <- for_permutation # whether used for perm or not
  
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
  #-----------------------------------------------------------------------------
  # score 1: proportions of 0
  
  # extract which nr is higher; nr of clusters or number of cell types  
  if(nr_clusters == 1){
    score_1 <- 0
  }else if(nr_clusters > 1){
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
  }
  
  print("score_1")
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # score 2: mean proportion of cells of a cell type per cluster
  
  mat_per_cluster[mat_per_cluster == 0] <- NA
  per_cluster_mean <- base::mean(mat_per_cluster, na.rm = TRUE)
  
  score_2 <- per_cluster_mean

  print("score_2")
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # score 3: mean purity of new clusters
  
  res_recl <- bluster::neighborPurity(mat_pca, 
                                      clusters = cluster_vector2, 
                                      k = 50)
  score_3 <- base::mean(res_recl$purity)
  
  print("score_3")
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # scores 4 - 6: established metrics for comparing two clusterings
  
  # Adjusted Rand index
  score_4 <- mclust::adjustedRandIndex(cluster_vector1, cluster_vector2)
  
  # Fowlkes-Mallows Index
  score_5 <- dendextend::FM_index_R(cluster_vector1, cluster_vector2)[1]
  
  # Variation of Information
  score_6 <- mcclust::vi.dist(S4Vectors::unfactor(cluster_vector1),
                              S4Vectors::unfactor(cluster_vector2))

  print("score_4 to score_6")
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  
  # depending on what this function is used for, output differs 
  if(!for_permutation){
    
    return_df <- base::data.frame(
      "type" = c("proportion_of_zeros", 
                 "mean_prop_cells_cluster",
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
      "mean_prop_cells_cluster" = score_2,
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

#################### permuting_reclustering_scores_seurat ######################

#-------------------------------------------------------------------------------

# function that combines two functions, one for re-clustering, one for 
# calculating scores

# used for permutation test with n iterations 
# used for seurat objects only
# the function will only output the scores for direct comparison, but not 
# the re-clustered datasets.

permuting_reclustering_scores_seurat <- function(
    iteration, 
    seu, 
    iteration_df, 
    resolution,
    data_use){
  
  print("starting iterations")
  
  iteration <- iteration # list of iterations like: as.list(c(1:i))
  seu <- seu # seurat object to be reclustered 
  iteration_df <- iteration_df # df with ncol=i and random gene positions/col
  resolution <- resolution # resolution to parse to FindClusters
  data_use <- data_use # data for reclustering ("log_counts" or "raw_counts")
  
  # subset object to previously generated random list of genes from iteration_df
  iteration_vector <- iteration_df[,iteration]

  seu_sub <- BiocGenerics::subset(
    seu,
    features = SeuratObject::Features(seu)[iteration_vector], 
    slot = "count")
  # check that GENES were subsetted
  print(nrow(seu_sub))

  print("starting standard_seu_pipeline function")
  # re-cluster seurat objects from using standard seurat pipeline
  seu_rec <- standard_seu_pipeline(
    seu = seu_sub, 
    features = SeuratObject::Features(seu_sub), 
    resolution = resolution,
    data_use = data_use,
    calc_umap = FALSE)

  print("starting calculate_scores function")
  # calculate re-clustering scores 
  res_df <- calculate_scores(seu_rec, for_permutation = TRUE)
  res_df$iteration[1] <- iteration
  
  return(res_df)
}



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

##################### permuting_reclustering_scores_sce ########################

#-------------------------------------------------------------------------------

# function that combines two functions, one for re-clustering, one for 
# calculating scores

# used for permutation test with n iterations 
# used for SCE objects (own datasets) only
# the function will only output the scores for direct comparison, but not 
# the re-clustered datasets.

permuting_reclustering_scores_sce <- function(
  iteration, 
  sce, 
  iteration_df, 
  resolution_louvain,
  k_graph){
  
  print("starting iteration")
  
  iteration <- iteration # list of iterations like: as.list(c(1:i))
  sce <- sce # sce object (own dataset) to be reclustered 
  iteration_df <- iteration_df # df with ncol=i and random gene positions/col 
  resolution_louvain <- resolution_louvain # resolution for louvain clustering
  k_graph <- k_graph # number of k for SNNgraph
  
  # subset object to previously generated random list of genes from iteration_df
  iteration_vector <- iteration_df[,iteration]

  # subset genes (rows)
  sce_sub <- sce[iteration_vector,]
  
  # make sure that only GENES have been subsetted, not cells
  print(nrow(sce))
  print(ncol(sce))
  
  print("starting reclustering function")
  # re-cluster sce objects using original pipeline
  seu_rec <- clustering_orig(
    sce = sce_sub,
    k_graph = k_graph,
    resolution_louvain = resolution_louvain
  )
  
  print("starting calculate_scores function")
  # calculate re-clustering scores 
  
  # define input
  cluster_vector1 <- seu_rec$celltypes
  cluster_vector2 <- seu_rec$reclustered
  mat_pca <- SingleCellExperiment::reducedDims(seu_rec)$PCA[,1:10]
  
  res_df <- calculate_scores( 
    cluster_vector1 = cluster_vector1,
    cluster_vector2 = cluster_vector2,
    mat_pca = mat_pca, 
    for_permutation = TRUE)
  
  return(res_df)
}

