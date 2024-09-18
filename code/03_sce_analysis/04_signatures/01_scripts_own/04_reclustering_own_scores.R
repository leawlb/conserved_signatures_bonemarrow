
# calculate scores for reclustered own datasets

#-------------------------------------------------------------------------------

set.seed(37)

library(mclust, quietly = TRUE)
library(mcclust, quietly = TRUE)
library(bluster, quietly = TRUE)
library(dendextend, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)

source(file = snakemake@params[["functions_reclustering"]])

#-------------------------------------------------------------------------------
# load SCE object with reclustered populations

sce <- base::readRDS(snakemake@input[["sce_input"]])

# cts to remove
cts_exclude <- snakemake@params[["cts_exclude"]]
fraction_curr <- snakemake@wildcards[["fraction"]]

# remove cts (failsave) to exclude and reorder factors for nicer plots
sce <- sce[,!sce$celltypes %in% cts_exclude]
sce$celltypes <- factor(
  sce$celltypes,
  levels = levels(sce$celltypes)[levels(sce$celltypes) %in% sce$celltypes])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# calculate scores

# signature 
cluster_vector1 <- sce$celltypes # cluster_vector1 = cell types
cluster_vector2 <- sce$cluster_signt # cluster_vector2 = new clusters
mat_pca <- SingleCellExperiment::reducedDims(sce)$PCA[,1:10]

# use own function for calculation
res_df_sign <- calculate_scores( 
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_sign$conservation_level <- base::rep(
  "conserved_signature", 
  nrow(res_df_sign))
res_df_sign$nr_celltypes <- base::rep(
  length(base::unique(cluster_vector1)), 
  nrow(res_df_sign))
res_df_sign$nr_clusters <- base::rep(
  length(base::unique(cluster_vector2)), 
  nrow(res_df_sign))
res_df_sign$fraction <- base::rep(
  fraction_curr, 
  nrow(res_df_sign))
res_df_sign$nr_genes_used <- base::rep(
  sce$cluster_signt_genes_used[1],
  nrow(res_df_sign))

#-------------------------------------------------------------------------------

# conserved markers 
cluster_vector2 <- sce$cluster_consm

res_df_mark <- calculate_scores(
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_mark$conservation_level <- base::rep("conserved_markers", 
                                            nrow(res_df_mark))

res_df_mark$nr_celltypes <- base::rep(
  length(base::unique(cluster_vector1)), 
  nrow(res_df_mark))
res_df_mark$nr_clusters <- base::rep(
  length(base::unique(cluster_vector2)), 
  nrow(res_df_mark))
res_df_mark$fraction <- base::rep(
  fraction_curr, 
  nrow(res_df_mark))
res_df_mark$nr_genes_used <- base::rep(
  sce$cluster_consm_genes_used[1],
  nrow(res_df_mark))

#-------------------------------------------------------------------------------

# all bl6 markers
cluster_vector2 <- sce$cluster_mmusm

res_df_mmsm <- calculate_scores(
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_mmsm$conservation_level <- base::rep(
  "mmusall_markers", 
  nrow(res_df_mmsm))
res_df_mmsm$nr_celltypes <- base::rep(
  length(base::unique(cluster_vector1)), 
  nrow(res_df_mmsm))
res_df_mmsm$nr_clusters <- base::rep(
  length(base::unique(cluster_vector2)), 
  nrow(res_df_mmsm))
res_df_mmsm$fraction <- base::rep(
  fraction_curr, 
  nrow(res_df_mmsm))
res_df_mmsm$nr_genes_used <- base::rep(
  sce$cluster_mmusm_genes_used[1],
  nrow(res_df_mmsm))

#-------------------------------------------------------------------------------

# ndges
cluster_vector2 <- sce$cluster_ndges

res_df_ndge <- calculate_scores(
  cluster_vector1 = cluster_vector1,
  cluster_vector2 = cluster_vector2,
  mat_pca = mat_pca)

res_df_ndge$conservation_level <- base::rep(
  "ndges", 
  nrow(res_df_ndge))
res_df_ndge$nr_celltypes <- base::rep(
  length(base::unique(cluster_vector1)), 
  nrow(res_df_ndge))
res_df_ndge$nr_clusters <- base::rep(
  length(base::unique(cluster_vector2)), 
  nrow(res_df_ndge))
res_df_ndge$fraction <- base::rep(
  fraction_curr,
  nrow(res_df_ndge))
res_df_ndge$nr_genes_used <- base::rep(
  sce$cluster_ndges_genes_used[1],
  nrow(res_df_ndge))

#-------------------------------------------------------------------------------

# combine
res_df_all <- base::rbind(res_df_sign, res_df_mark, res_df_mmsm, res_df_ndge)

#-------------------------------------------------------------------------------
base::saveRDS(res_df_all, snakemake@output[["score_df"]])

utils::sessionInfo()
