
#-------------------------------------------------------------------------------

# prepare datasets from tabula sapiens for downstream analysis
# three datasets from TS: bone marrow, hspc_progenitors, and stromal cells
# bone marrow (immune cells) from adult donors
# hspc_progenitors and stromal cells from embryonal donors

set.seed(37)

library(Seurat, quietly = TRUE)

#-------------------------------------------------------------------------------
# get paths from snakemake

ts_hsc_progenitors_input <- snakemake@input[["ts_hsc_progenitors_input"]]
ts_bone_marrow_input <- snakemake@input[["ts_bone_marrow_input"]]
ts_stromal_input <- snakemake@input[["ts_stromal_input"]]

ts_hsc_progenitors_output<- snakemake@output[["ts_hscs_progenitors"]]
ts_bone_marrow_output <- snakemake@output[["ts_bone_marrow"]]
ts_stromal_output <- snakemake@output[["ts_all_stromal"]]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# HSCs/PROGENITORS

seu_ts_hscs_progenitors <- base::readRDS(ts_hsc_progenitors_input)

#-------------------------------------------------------------------------------
# check sex and age of donors (+ other data)
# these are embryonic hematopoietic cells from Carnegie stage 13 to
# 17th week post-fertilization human stage
# some of the cells were derived from yolk sac etc.

print(base::table(seu_ts_hscs_progenitors$sex))
print(base::table(seu_ts_hscs_progenitors$tissue))
print(base::table(seu_ts_hscs_progenitors$development_stage))
print(base::table(seu_ts_hscs_progenitors$is_maternal_contaminant))
print(base::table(seu_ts_hscs_progenitors$disease))
print(base::table(seu_ts_hscs_progenitors$organism))
print(base::table(seu_ts_hscs_progenitors$assay))

# print(base::table(seu_ts_hscs_progenitors$celltype_annotation,
#                   seu_ts_hscs_progenitors$cell_type))

#-------------------------------------------------------------------------------

# remove cell types not in our dataset, as we cannot separate them
cts_to_remove <- c(
  "fraction A pre-pro B cell",
  "double negative thymocyte",
  "megakaryocyte",
  "unknown"
)

seu_ts_hscs_progenitors <- seu_ts_hscs_progenitors[,!seu_ts_hscs_progenitors$cell_type %in% cts_to_remove]

# re-factor cell types to desired order
# keep in "cell_type" slot for downstream compatibility
seu_ts_hscs_progenitors$cell_type <- factor(
  seu_ts_hscs_progenitors$cell_type, 
  levels = c(
    "hematopoietic stem cell",
    "hematopoietic multipotent progenitor cell",
    "common myeloid progenitor",
    "granulocyte monocyte progenitor cell",
    "promonocyte",
    "promyelocyte",
    "early lymphoid progenitor",
    "megakaryocyte-erythroid progenitor cell"
    ))

print(base::table(is.na(seu_ts_hscs_progenitors$cell_type)))
stopifnot(!is.na(seu_ts_hscs_progenitors$cell_type))

# print(base::table(seu_ts_hscs_progenitors$celltype_annotation,
#                   seu_ts_hscs_progenitors$cell_type))

#-------------------------------------------------------------------------------

# remove umap coordinates and neighbors
print(seu_ts_hscs_progenitors@reductions)
# no PCA
seu_ts_hscs_progenitors@neighbors <- list()
seu_ts_hscs_progenitors@reductions <- list()

dim(seu_ts_hscs_progenitors)

#-------------------------------------------------------------------------------

# add info on which column of the ensembl data frame to use based on Features 
seu_ts_hscs_progenitors@misc$ensembl_column_use <- "ENSG_ID" # human IDs

#-------------------------------------------------------------------------------
# add info on which assay should be used for reclustering
seu_ts_hscs_progenitors@misc$data_use <- "raw_counts"

#-------------------------------------------------------------------------------

base::saveRDS(seu_ts_hscs_progenitors, ts_hsc_progenitors_output)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# BONE MARROW

seu_ts_bone_marrow <- base::readRDS(ts_bone_marrow_input)

#-------------------------------------------------------------------------------
# check sex, age, and other data

print(base::table(seu_ts_bone_marrow$sex))
print(base::table(seu_ts_bone_marrow$development_stage))
print(base::table(seu_ts_bone_marrow$disease))
print(base::table(seu_ts_bone_marrow$organism))
print(base::table(seu_ts_bone_marrow$assay))

# print(base::table(seu_ts_bone_marrow$free_annotation,
#                   seu_ts_bone_marrow$cell_type))

#-------------------------------------------------------------------------------

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
# keep in "cell_type" slot for downstream compatibility
seu_ts_bone_marrow$cell_type <- factor(seu_ts_bone_marrow$cell_type,
                                        levels = c(
                                          "hematopoietic stem cell",
                                          "common myeloid progenitor",
                                          "erythroid progenitor cell",
                                          "granulocyte"
                                        ))

stopifnot(!is.na(seu_ts_bone_marrow$cell_type))

# print(base::table(seu_ts_bone_marrow$free_annotation,
#                   seu_ts_bone_marrow$cell_type))

#-------------------------------------------------------------------------------

# remove umap coordinates and neighbors
seu_ts_bone_marrow@neighbors <- list()
# keep original PCA coordinates for comparison later but remove other reductions
seu_ts_bone_marrow_pca <- seu_ts_bone_marrow
seu_ts_bone_marrow@reductions <- list()
seu_ts_bone_marrow@reductions$pca_orig <- seu_ts_bone_marrow_pca@reductions$pca

print("after removal")
print(seu_ts_bone_marrow@reductions)
dim(seu_ts_bone_marrow)

#-------------------------------------------------------------------------------

# add info on which column of the ensembl data frame to use based on Features 
seu_ts_bone_marrow@misc$ensembl_column_use <- "ENSG_ID" # human IDs

#-------------------------------------------------------------------------------
# add info on which assay should be used for reclustering
seu_ts_bone_marrow@misc$data_use <- "raw_counts"

#-------------------------------------------------------------------------------

base::saveRDS(seu_ts_bone_marrow, ts_bone_marrow_output)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# STROMAL 

ts_all_stromal <- base::readRDS(ts_stromal_input)

#-------------------------------------------------------------------------------
# check sex, age, and other data
# also embryonal

print(base::table(ts_all_stromal$sex))
print(base::table(ts_all_stromal$development_stage))
print(base::table(ts_all_stromal$disease))
print(base::table(ts_all_stromal$organism))
print(base::table(ts_all_stromal$assay))
print(base::table(ts_all_stromal$is_maternal_contaminant))
print(base::table(ts_all_stromal$tissue))

# print(base::table(ts_all_stromal$celltype_annotation,
#                   ts_all_stromal$cell_type))

#-------------------------------------------------------------------------------

# keep only cells from bone marrow
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
# keep in "cell_type" slot
ts_all_stromal$cell_type <- factor(ts_all_stromal$cell_type,
                                   levels = c(
                                     "osteoblast",
                                     "chondrocyte",
                                     "fibroblast",
                                     "myofibroblast cell",
                                     "endothelial cell",
                                     "vascular associated smooth muscle cell"
                                   ))

stopifnot(!is.na(ts_all_stromal$cell_type))

# print(base::table(ts_all_stromal$celltype_annotation,
#                   ts_all_stromal$cell_type))

#-------------------------------------------------------------------------------

# remove umap coordinates and neighbors
print(ts_all_stromal@reductions)
# no PCA
ts_all_stromal@neighbors <- list()
ts_all_stromal@reductions <- list()

dim(ts_all_stromal)

#-------------------------------------------------------------------------------

# add info on which column of the ensembl data frame to use based on Features 
ts_all_stromal@misc$ensembl_column_use <- "ENSG_ID" # human IDs

#-------------------------------------------------------------------------------
# add info on which assay should be used for reclustering
ts_all_stromal@misc$data_use <- "raw_counts"

#-------------------------------------------------------------------------------

base::saveRDS(ts_all_stromal, ts_stromal_output)

utils::sessionInfo()