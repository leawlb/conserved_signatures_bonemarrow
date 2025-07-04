#-------------------------------------------------------------------------------

library(scran, quietly = TRUE)
library(scmap, quietly = TRUE)
set.seed(37)

# From https://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html

sce <- base::readRDS(file = snakemake@input[["sce_input"]])

# load reference datasets 
# reference datasets have been downloaded, pre-processed, and stored in
# metadata/scRNAseq/01_sce_prep/references_raw and .../references_processed
ref_baccin <- base::readRDS(file = snakemake@input[["ref_baccin_sce"]])
ref_dahlin <- base::readRDS(file = snakemake@input[["ref_dahlin_sce"]])
ref_dolgalev <- base::readRDS(file = snakemake@input[["ref_dolgalev_sce"]])

#-------------------------------------------------------------------------------

# Baccin et al.

# prepare reference and test datasets
rowData(ref_baccin)$feature_symbol <- rownames(ref_baccin)
# get most informative features for projection
ref_baccin <- scmap::selectFeatures(ref_baccin, suppress_plot = TRUE)

test_sce <- sce
rowData(test_sce)$feature_symbol <- rownames(test_sce)

# get index lists (required for projection)
# calculates centroids of each cell type and merge them into a single table
# for cluster projection
ref_baccin_clst_ind <- scmap::indexCluster(ref_baccin, 
                                           cluster_col = "identity_ref")
# enables fast approximate nearest neighbor search for cell projection
ref_baccin_cell_ind <- scmap::indexCell(ref_baccin)

# scmapcell
# for each cell in test dataset, search the nearest neighbor in ref dataset
set.seed(37)
ref_baccin_cell_temp <- scmap::scmapCell(
  projection = test_sce, 
  index_list = list(results = metadata(ref_baccin_cell_ind)$scmap_cell_index))

# cell type classification from scmapcell results/projection
ref_baccin_cell_results <- scmap::scmapCell2Cluster(
  scmapCell_results = ref_baccin_cell_temp, 
  cluster_list = list(as.character(colData(ref_baccin)$identity_ref)))

# scmapclust
# projects cells from test dataset into ref dataset and obtains labels
ref_baccin_clst_results <- scmap::scmapCluster(
  projection = test_sce,
  index_list = list(metadata(ref_baccin_clst_ind)$scmap_cluster_index), 
  threshold = 0.7)

sce$baccin_celltype_scmapcell <- ref_baccin_cell_results$combined_labs
sce$baccin_celltype_scmapcell_sim <- ref_baccin_cell_results$scmap_cluster_siml

sce$baccin_celltype_scmapclust <- ref_baccin_clst_results$combined_labs
sce$baccin_celltype_scmapclust_sim <- ref_baccin_clst_results$scmap_cluster_siml

#-------------------------------------------------------------------------------

# Dahlin et al.

rowData(ref_dahlin)$feature_symbol <- rownames(ref_dahlin)
ref_dahlin <- scmap::selectFeatures(ref_dahlin, suppress_plot = TRUE)

test_sce <- sce
rowData(test_sce)$feature_symbol <- rownames(test_sce)

ref_dahlin_clst_ind <- scmap::indexCluster(ref_dahlin,
                                           cluster_col = "identity_ref")
ref_dahlin_cell_ind <- scmap::indexCell(ref_dahlin)

set.seed(37)
ref_dahlin_cell_temp <- scmap::scmapCell(
  projection = test_sce, 
  index_list = list(results = metadata(ref_dahlin_cell_ind)$scmap_cell_index))

ref_dahlin_cell_results <- scmap::scmapCell2Cluster(
  scmapCell_results = ref_dahlin_cell_temp, 
  cluster_list = list(as.character(colData(ref_dahlin)$identity_ref)))

ref_dahlin_clst_results <- scmap::scmapCluster(
  projection = test_sce,
  index_list = list(metadata(ref_dahlin_clst_ind)$scmap_cluster_index), 
  threshold = 0.7)

sce$dahlin_celltype_scmapcell <- ref_dahlin_cell_results$combined_labs
sce$dahlin_celltype_scmapcell_sim <- ref_dahlin_cell_results$scmap_cluster_siml

sce$dahlin_celltype_scmapclust <- ref_dahlin_clst_results$combined_labs
sce$dahlin_celltype_scmapclust_sim <- ref_dahlin_clst_results$scmap_cluster_siml

#-------------------------------------------------------------------------------

# Dolgalev et al.

rowData(ref_dolgalev)$feature_symbol <- rownames(ref_dolgalev)
ref_dolgalev <- scmap::selectFeatures(ref_dolgalev, suppress_plot = TRUE)

test_sce <- sce
rowData(test_sce)$feature_symbol <- rownames(test_sce)

ref_dolgalev_clst_ind <- scmap::indexCluster(ref_dolgalev, 
                                             cluster_col = "labelsimple")
ref_dolgalev_cell_ind <- scmap::indexCell(ref_dolgalev)

set.seed(37)
ref_dolgalev_cell_temp <- scmap::scmapCell(
  projection = test_sce, 
  index_list = list(results = metadata(ref_dolgalev_cell_ind)$scmap_cell_index))

ref_dolgalev_cell_results <- scmap::scmapCell2Cluster(
  scmapCell_results = ref_dolgalev_cell_temp, 
  cluster_list = list(as.character(colData(ref_dolgalev)$labelsimple)))

ref_dolgalev_clst_results <- scmap::scmapCluster(
  projection = test_sce,
  index_list = list(metadata(ref_dolgalev_clst_ind)$scmap_cluster_index), 
  threshold = 0.7)

sce$dolgalev_celltype_scmapcell <- ref_dolgalev_cell_results$combined_labs
sce$dolgalev_celltype_scmapcell_sim <- ref_dolgalev_cell_results$scmap_cluster_siml

sce$dolgalev_celltype_scmapclust <- ref_dolgalev_clst_results$combined_labs
sce$dolgalev_celltype_scmapclust_sim <- ref_dolgalev_clst_results$scmap_cluster_siml

#-------------------------------------------------------------------------------

print(sce)
base::saveRDS(sce, file = snakemake@output[["sce_output"]])

utils::sessionInfo()