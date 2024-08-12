
#-------------------------------------------------------------------------------

# prepare datasets from tabula sapiens for downstream analysis

set.seed(37)

library(Seurat, quietly = TRUE)

#-------------------------------------------------------------------------------
# get paths from snakemake

ts_bone_marrow_input <- snakemake@input[["ts_bone_marrow_input"]]
ts_hsc_progenitors_input <- snakemake@input[["ts_hsc_progenitors_input"]]
ts_stromal_input <- snakemake@input[["ts_stromal_input"]]

ts_bone_marrow_output <- snakemake@output[["ts_bone_marrow_output"]]
ts_hsc_progenitors_output<- snakemake@output[["ts_hsc_progenitors_output"]]
ts_stromal_output <- snakemake@output[["ts_stromal_output"]]

#-------------------------------------------------------------------------------

seu_ts_hscs_progenitors <- base::readRDS(ts_bone_marrow_input)

cts_to_remove <- c(
  "fraction A pre-pro B cell",
  "double negative thymocyte",
  "megakaryocyte",
  "unknown"
)

seu_ts_hscs_progenitors <- seu_ts_hscs_progenitors[,!seu_ts_hscs_progenitors$cell_type %in% cts_to_remove]

# re-factor cell types to desired order
# will also keep cells that are not from tissue "bone marrow"
seu_ts_hscs_progenitors$cell_type <- factor(seu_ts_hscs_progenitors$cell_type, 
                                            levels = c(
                                              "hematopoietic stem cell",
                                              "hematopoietic multipotent progenitor cell",
                                              "common myeloid progenitor",
                                              "granulocyte monocyte progenitor cell",
                                              "promyelocyte",
                                              "promonocyte",
                                              "early lymphoid progenitor",
                                              "megakaryocyte-erythroid progenitor cell"
                                            ))

# remove umap coordinates and neighbors
seu_ts_hscs_progenitors@neighbors <- list()
# keep original PCA coordinates for comparison later but remove other reductions
seu_ts_hscs_progenitors_pca <- seu_ts_hscs_progenitors
seu_ts_hscs_progenitors@reductions <- list()
seu_ts_hscs_progenitors@reductions$pca_orig <- seu_ts_hscs_progenitors_pca@reductions$pca

print("after removal")
print(seu_ts_hscs_progenitors@reductions)
dim(seu_ts_hscs_progenitors)

# add info on which column of the ensembl data frame to use based on Features 
seu_ts_hscs_progenitors@misc$ensembl_column_use <- "ENSG_ID" # human IDs

base::saveRDS(seu_ts_hscs_progenitors, ts_hsc_progenitors_output)

#-------------------------------------------------------------------------------

seu_ts_bone_marrow <- base::readRDS(ts_bone_marrow_input)

# remove mature cell types
cts_to_remove <- c(
  "CD4-positive, alpha-beta T cell",
  "CD8-positive, alpha-beta T cell",
  "mature NK T cell",
  "erythrocyte",
  "macrophage",
  "naive B cell",
  "memory B cell",
  "plasma cell",
  "plasmablast",
  "neutrophil",
  "monocyte"
)

seu_ts_bone_marrow <- seu_ts_bone_marrow[,!seu_ts_bone_marrow$cell_type %in% cts_to_remove]

# re-factor cell_type
seu_ts_bone_marrow$cell_type <- factor(seu_ts_bone_marrow$cell_type,
                                        levels = c(
                                          "hematopoietic stem cell",
                                          "common myeloid progenitor",
                                          "granulocyte",
                                          "erythroid progenitor cell"
                                        ))

# remove umap coordinates and neighbors
seu_ts_bone_marrow@neighbors <- list()
# keep original PCA coordinates for comparison later but remove other reductions
seu_ts_bone_marrow_pca <- seu_ts_bone_marrow
seu_ts_bone_marrow@reductions <- list()
seu_ts_bone_marrow@reductions$pca_orig <- seu_ts_bone_marrow_pca@reductions$pca

print("after removal")
print(seu_ts_bone_marrow@reductions)

dim(seu_ts_bone_marrow)

# add info on which column of the ensembl data frame to use based on Features 
seu_ts_bone_marrow@misc$ensembl_column_use <- "ENSG_ID" # human IDs

base::saveRDS(seu_ts_bone_marrow, ts_bone_marrow_output)

#-------------------------------------------------------------------------------

ts_all_stromal <- base::readRDS(ts_stromal_input)

# keep only cell originated from bone marrow
ts_all_stromal <- ts_all_stromal[,ts_all_stromal$tissue == "bone marrow"]

# remove the cell types with only very small nr of cells left
ts_all_stromal <- ts_all_stromal[,!ts_all_stromal$cell_type %in% c(
  "epithelial cell of nephron",
  "hepatocyte",
  "eurydendroid cell",
  "unknown",
  "cell of skeletal muscle",
  "skeletal muscle satellite cell"
)]

# re-factor cell_type
ts_all_stromal$cell_type <- factor(ts_all_stromal$cell_type,
                                   levels = c(
                                     "osteoblast",
                                     "chondrocyte",
                                     "fibroblast",
                                     "myofibroblast cell",
                                     "endothelial cell",
                                     "vascular associated smooth muscle cell"
                                   ))

# remove umap coordinates and neighbors
ts_all_stromal@neighbors <- list()
# keep original PCA coordinates for comparison later but remove other reductions
ts_all_stromal_pca <- ts_all_stromal
print("check stromal before removal")
print(ts_all_stromal_pca@reductions)
ts_all_stromal@reductions <- list()
ts_all_stromal@reductions$pca_orig <- ts_all_stromal_pca@reductions$pca

print("after removal")
print(ts_all_stromal@reductions)

dim(ts_all_stromal)

# add info on which column of the ensembl data frame to use based on Features 
ts_all_stromal@misc$ensembl_column_use <- "ENSG_ID" # human IDs

base::saveRDS(ts_all_stromal, ts_stromal_output)

utils::sessionInfo()