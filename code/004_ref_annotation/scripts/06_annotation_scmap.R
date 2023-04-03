#-------------------------------------------------------------------------------

library(DropletUtils)
library(scmap)
set.seed(37)

# From https://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html

sce <- readRDS(file = snakemake@input[["sce_05"]])

# load reference datasets 
ref_baccin <- readRDS(file = snakemake@input[["ref_baccin_sce"]])
ref_dahlin <- readRDS(file = snakemake@input[["ref_dahlin_sce"]])
ref_dolgalev <- readRDS(file = snakemake@input[["ref_dolgalev_sce"]])

#-------------------------------------------------------------------------------

# Baccin

# prepare reference and test datasets
rowData(ref_baccin)$feature_symbol <- rownames(ref_baccin)
ref_baccin <- selectFeatures(ref_baccin, suppress_plot = TRUE)

test_sce <- sce
rowData(test_sce)$feature_symbol <- rownames(test_sce)

# get index lists
ref_baccin_clst_ind <- indexCluster(ref_baccin, cluster_col = "identity_ref")
ref_baccin_cell_ind <- indexCell(ref_baccin)

# scmapcell
set.seed(37)
ref_baccin_cell_temp <- scmapCell(
  projection = test_sce, 
  index_list = list(results = metadata(ref_baccin_cell_ind)$scmap_cell_index))

ref_baccin_cell_results <- scmapCell2Cluster(
  scmapCell_results = ref_baccin_cell_temp, 
  cluster_list = list(as.character(colData(ref_baccin)$identity_ref)))

# scmapclust
ref_baccin_clst_results <- scmapCluster(
  projection = test_sce,
  index_list = list(metadata(ref_baccin_clst_ind)$scmap_cluster_index), 
  threshold = 0.7)

sce$baccin_celltype_scmapcell <- ref_baccin_cell_results$combined_labs
sce$baccin_celltype_scmapcell_sim <- ref_baccin_cell_results$scmap_cluster_siml

sce$baccin_celltype_scmapclust <- ref_baccin_clst_results$combined_labs
sce$baccin_celltype_scmapclust_sim <- ref_baccin_clst_results$scmap_cluster_siml

#-------------------------------------------------------------------------------

# Dahlin

# prepare reference and test datasets
rowData(ref_dahlin)$feature_symbol <- rownames(ref_dahlin)
ref_dahlin <- selectFeatures(ref_dahlin, suppress_plot = TRUE)

test_sce <- sce
rowData(test_sce)$feature_symbol <- rownames(test_sce)

# get index lists
ref_dahlin_clst_ind <- indexCluster(ref_dahlin, cluster_col = "identity_ref")
ref_dahlin_cell_ind <- indexCell(ref_dahlin)

# scmapcell
set.seed(37)
ref_dahlin_cell_temp <- scmapCell(
  projection = test_sce, 
  index_list = list(results = metadata(ref_dahlin_cell_ind)$scmap_cell_index))

ref_dahlin_cell_results <- scmapCell2Cluster(
  scmapCell_results = ref_dahlin_cell_temp, 
  cluster_list = list(as.character(colData(ref_dahlin)$identity_ref)))

# scmapclust
ref_dahlin_clst_results <- scmapCluster(
  projection = test_sce,
  index_list = list(metadata(ref_dahlin_clst_ind)$scmap_cluster_index), 
  threshold = 0.7)

sce$dahlin_celltype_scmapcell <- ref_dahlin_cell_results$combined_labs
sce$dahlin_celltype_scmapcell_sim <- ref_dahlin_cell_results$scmap_cluster_siml

sce$dahlin_celltype_scmapclust <- ref_dahlin_clst_results$combined_labs
sce$dahlin_celltype_scmapclust_sim <- ref_dahlin_clst_results$scmap_cluster_siml

#-------------------------------------------------------------------------------

# Baccin

# prepare reference and test datasets
rowData(ref_dolgalev)$feature_symbol <- rownames(ref_dolgalev)
ref_dolgalev <- selectFeatures(ref_dolgalev, suppress_plot = TRUE)

test_sce <- sce
rowData(test_sce)$feature_symbol <- rownames(test_sce)

# get index lists
ref_dolgalev_clst_ind <- indexCluster(ref_dolgalev, cluster_col = "labelsimple")
ref_dolgalev_cell_ind <- indexCell(ref_dolgalev)

# scmapcell
set.seed(37)
ref_dolgalev_cell_temp <- scmapCell(
  projection = test_sce, 
  index_list = list(results = metadata(ref_dolgalev_cell_ind)$scmap_cell_index))

ref_dolgalev_cell_results <- scmapCell2Cluster(
  scmapCell_results = ref_dolgalev_cell_temp, 
  cluster_list = list(as.character(colData(ref_dolgalev)$labelsimple)))

# scmapclust
ref_dolgalev_clst_results <- scmapCluster(
  projection = test_sce,
  index_list = list(metadata(ref_dolgalev_clst_ind)$scmap_cluster_index), 
  threshold = 0.7)

sce$dolgalev_celltype_scmapcell <- ref_dolgalev_cell_results$combined_labs
sce$dolgalev_celltype_scmapcell_sim <- ref_dolgalev_cell_results$scmap_cluster_siml

sce$dolgalev_celltype_scmapclust <- ref_dolgalev_clst_results$combined_labs
sce$dolgalev_celltype_scmapclust_sim <- ref_dolgalev_clst_results$scmap_cluster_siml

print(sce)

#-------------------------------------------------------------------------------

# saves individuals in species folders
saveRDS(sce, file = snakemake@output[["sce_06"]])
