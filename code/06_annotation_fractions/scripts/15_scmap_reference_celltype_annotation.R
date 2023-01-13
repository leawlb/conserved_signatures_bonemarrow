#-------------------------------------------------------------------------------

library(DropletUtils)
library(scmap)
set.seed(37)

# From https://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html

sce <- readRDS(file = snakemake@input[["sce_14"]])

ref_baccin <- readRDS(file = snakemake@params[["ref_baccin_sce"]])
ref_dahlin <- readRDS(file = snakemake@params[["ref_dahlin_sce"]])
ref_dolgalev <- readRDS(file = snakemake@params[["ref_dolgalev_sce"]])

# Baccin
rowData(ref_baccin)$feature_symbol <- rownames(ref_baccin)
ref_baccin <- selectFeatures(ref_baccin, suppress_plot = TRUE)
ref_baccin <- indexCluster(ref_baccin, cluster_col = "identity_ref")

temp_sce <- sce
rowData(temp_sce)$feature_symbol <- rownames(temp_sce)
results <- scmapCluster(projection = temp_sce,
                        index_list = list(
                          metadata(ref_baccin)$scmap_cluster_index), 
                        threshold = 0.7)

sce$baccin_identity_pred_cluster <- results$combined_labs
sce$baccin_identity_pred_cluster_sim <- results$scmap_cluster_siml

# do not separate by fraction anymore because scmap can reject
# Dahlin
rowData(ref_dahlin)$feature_symbol <- rownames(ref_dahlin)
ref_dahlin <- selectFeatures(ref_dahlin, suppress_plot = TRUE)
ref_dahlin <- indexCluster(ref_dahlin, cluster_col = "identity_ref")

temp_sce <- sce
rowData(temp_sce)$feature_symbol <- rownames(temp_sce)
results <- scmapCluster(projection = temp_sce,
                        index_list = list(
                          metadata(ref_dahlin)$scmap_cluster_index),
                        threshold = 0.7)

sce$dahlin_identity_pred_cluster <- results$combined_labs
sce$dahlin_identity_pred_cluster_sim <- results$scmap_cluster_siml

# Dolgalev
rowData(ref_dolgalev)$feature_symbol <- rownames(ref_dolgalev)
ref_dolgalev <- selectFeatures(ref_dolgalev, suppress_plot = TRUE)
ref_dolgalev <- indexCluster(ref_dolgalev, cluster_col = "labelsimple")

temp_sce <- sce
rowData(temp_sce)$feature_symbol <- rownames(temp_sce)
results <- scmapCluster(projection = temp_sce,
                        index_list = list(
                          metadata(ref_dolgalev)$scmap_cluster_index), 
                        threshold = 0.7)

sce$dolgalev_identity_pred_cluster <- results$combined_labs
sce$dolgalev_identity_pred_cluster_sim <- results$scmap_cluster_siml

print(sce)
saveRDS(sce, file = snakemake@output[["sce_15"]])
