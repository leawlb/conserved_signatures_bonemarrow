#-------------------------------------------------------------------------------

# prepare datasets containing hematopoietic cells from naked mole rat (NMR)
# and zebrafish (ZEB)

set.seed(37)

library(Seurat, quietly = TRUE)
library(SeuratObject, quietly = TRUE)
library(stringr, quietly = TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ZEBRAFISH

# construct a Seurat object from the downloaded data containing all necessary
# assays and clustering information

#-------------------------------------------------------------------------------
# get all paths to the different files per assay

assays <- snakemake@input[["assays"]]
print(assays)

afnc <- "aggregated_filtered_normalised_counts"
afcs <- "aggregated_filtered_counts"
etpm <- "expression_tpm"

print(afnc)
print(afcs)
print(etpm)

#-------------------------------------------------------------------------------
# DATA = normalised counts

afnc_paths <- assays[base::grep(afnc, assays)]

counts_normalised <- Seurat::ReadMtx(
  mtx = afnc_paths[which(stringr::str_sub(afnc_paths, -4, -1) == ".mtx")],
  cells = afnc_paths[which(stringr::str_sub(afnc_paths, -4, -1) == "cols")],
  features = afnc_paths[which(stringr::str_sub(afnc_paths, -4, -1) == "rows")]
)
counts_normalised[1:4, 1:4]

#-------------------------------------------------------------------------------
# COUNTS = raw files 

afcs_paths <- assays[base::grep(afcs, assays)]

counts_raw <- Seurat::ReadMtx(
  mtx = afcs_paths[which(stringr::str_sub(afcs_paths, -4, -1) == ".mtx")],
  cells = afcs_paths[which(stringr::str_sub(afcs_paths, -4, -1) == "cols")],
  features = afcs_paths[which(stringr::str_sub(afcs_paths, -4, -1) == "rows")]
)
counts_raw[1:4, 1:4]

#-------------------------------------------------------------------------------
# METADATA = clusters

zeb_clustering_file <- snakemake@input[["zeb_clustering_file"]]

clusters <- utils::read.table(
  file = zeb_clustering_file,
  sep = '\t', 
  header = TRUE)
print(clusters[1:4,1:4])

# look at clustering info
selected_k <- clusters$K[clusters$sel.K == TRUE]
clusters <- clusters[,colnames(clusters) != "sel.K"]
clusters <- tibble::column_to_rownames(clusters, var = "K")
print(clusters[1:4,1:4])

rownames(clusters) <- base::paste0("K_", rownames(clusters))
metadata <- base::as.data.frame(t(clusters[,which(
  colnames(clusters) %in% colnames(counts_normalised))]))

metadata$selected_clustering <- metadata$K_10
print(head(metadata))

print(dim(metadata))
print(dim(counts_raw))
print(dim(counts_normalised))
#-------------------------------------------------------------------------------

# construct a seurat object from these files

seu_zeb <- SeuratObject::CreateSeuratObject(
  meta.data = metadata,
  counts = counts_raw, 
  data = counts_normalised,
  project = "SeuratProject", 
  assay = "RNA")

print(seu_zeb)
print("conversion done")

#-------------------------------------------------------------------------------
# add annotation

# I did not find the original cell type annotation in the available data
# (For now) I am annotating cell types ROUGHLY, based on gene expression
# of well-known marker genes

seu_zeb$manually_annotated_ct <- vector(length = ncol(seu_zeb))

seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 1] <- "Thrombo"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 2] <- "Neutro1"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 5] <- "Neutro2"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 3] <- "Active?"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 4] <- "Mono/Macro1"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 6] <- "Ery1"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 10] <- "Ery2"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 7] <- "Ery precursor"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 8] <- "Mono/Macro2"
seu_zeb$manually_annotated_ct[seu_zeb$selected_clustering == 9] <- "Lymphoid?"

print(seu_zeb)
print(head(seu_zeb@meta.data))

#-------------------------------------------------------------------------------

base::saveRDS(seu_zeb, snakemake@output[["zeb_seu"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# NAKED MOLE RAT

Seurat_hgl_sorted_BM <- snakemake@input[["Seurat_hgl_sorted_BM"]]
Seurat_hgl_whole_BM <- snakemake@input[["Seurat_hgl_whole_BM"]]

# load objects into environment, then rename
load(Seurat_hgl_sorted_BM)
seu_nmr_srt <- S.NOmp

load(Seurat_hgl_whole_BM)
seu_nmr_whl <- S.NOmp

# remove cell types that are not in our dataset/conserved signatures
seu_nmr_srt <- seu_nmr_srt[,!seu_nmr_srt$celltype.combi %in% 
                             c("TC", "PC", "ERY", "MO", "DC", "PB-GC", "BC")]

# will use sorted HSPCs for re-clustering, but also export whole HSPCs for
# testing and other uses
base::saveRDS(seu_nmr_srt, snakemake@output[["nmr_sorted_hspc"]])
base::saveRDS(seu_nmr_whl, snakemake@output[["nmr_whole_hspc"]])

utils::sessionInfo()
