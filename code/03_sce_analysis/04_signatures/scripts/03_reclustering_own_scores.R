
# calculate scores for reclustered own datasets

#-------------------------------------------------------------------------------

set.seed(37)

library(mclust, quietly = TRUE)
library(mcclust, quietly = TRUE)
library(bluster, quietly = TRUE)
library(cluster, quietly = TRUE)
library(dendextend, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)

source(file = snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
# load SCE object with reclustered populations

sce <- base::readRDS(snakemake@input[["sce_input"]])

cts_exclude <- snakemake@params[["cts_exclude"]]

# remove cts to exclude and reorder for nicer plots
sce <- sce[,!sce$celltypes %in% cts_exclude]
sce$celltypes <- factor(
  sce$celltypes,
  levels = levels(sce$celltypes)[levels(sce$celltypes) %in% sce$celltypes])

#-------------------------------------------------------------------------------
# calculate scores

# signature 
cluster_vector1 <- sce$celltypes
cluster_vector2 <- sce$cluster_signt
mat_pca <- SingleCellExperiment::reducedDims(sce)$PCA[,1:10]

res_df_sign <- calculate_scores_long( # own function
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_sign$conservation_level <- base::rep("conserved_signature", 
                                            nrow(res_df_sign))
    
# conserved markers 
cluster_vector2 <- sce$cluster_consm

res_df_mark <- calculate_scores_long(
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_mark$conservation_level <- base::rep("conserved_markers", 
                                            nrow(res_df_mark))
res_df_mark  

# all bl6 markers
cluster_vector2 <- sce$cluster_mmusm

res_df_mmsm <- calculate_scores_long(
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_mmsm$conservation_level <- base::rep("mmusall_markers", 
                                            nrow(res_df_mmsm))

# ndges
cluster_vector2 <- sce$cluster_ndges

res_df_ndge <- calculate_scores_long(
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_ndge$conservation_level <- base::rep("ndges", 
                                            nrow(res_df_ndge))
#-------------------------------------------------------------------------------

# combine
res_df_all <- base::rbind(res_df_sign, res_df_mark, res_df_mmsm, res_df_ndge)

#-------------------------------------------------------------------------------
base::saveRDS(res_df_all, snakemake@output[["score_df"]])

utils::sessionInfo()
