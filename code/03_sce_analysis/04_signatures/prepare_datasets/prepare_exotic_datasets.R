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
# this is how the unzipped files are called

assay_path <- snakemake@input[["assays"]]
print(assay_path)

afnc <- "aggregated_filtered_normalised_counts"
afcs <- "aggregated_filtered_counts"
etpm <- "expression_tpm"

#-------------------------------------------------------------------------------
# DATA = normalised counts

afnc_paths <- assay_path[base::grep(afnc, assay_path)]

counts_normalised <- Seurat::ReadMtx(
  mtx = afnc_paths[which(stringr::str_sub(afnc_paths, -4, -1) == ".mtx")],
  cells = afnc_paths[which(stringr::str_sub(afnc_paths, -4, -1) == "cols")],
  features = afnc_paths[which(stringr::str_sub(afnc_paths, -4, -1) == "rows")]
)
print(counts_normalised[1:4, 1:4])

#-------------------------------------------------------------------------------
# COUNTS = raw files 

afcs_paths <- assay_path[base::grep(afcs, assay_path)]

counts_raw <- Seurat::ReadMtx(
  mtx = afcs_paths[which(stringr::str_sub(afcs_paths, -4, -1) == ".mtx")],
  cells = afcs_paths[which(stringr::str_sub(afcs_paths, -4, -1) == "cols")],
  features = afcs_paths[which(stringr::str_sub(afcs_paths, -4, -1) == "rows")]
)
print(counts_raw[1:4, 1:4])

#-------------------------------------------------------------------------------
# METADATA = clusters

zeb_clustering_file <- snakemake@input[["zeb_clustering_file"]]

clusters <- utils::read.table(
  file = zeb_clustering_file,
  sep = '\t', 
  header = TRUE)
print(clusters[1:4,1:4])

# look at clustering info
# data was clustered using different ks, k = 10 is selected
# no other document with cell type annotation was found
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
# no reduced dimensions

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

# factorize cell types and put in correct slot for downstream analysis
seu_zeb$manually_annotated_ct <- factor(
  seu_zeb$manually_annotated_ct,
  levels = c("Thrombo",
             "Ery precursor",
             "Ery1",
             "Ery2",
             "Active?",
             "Neutro1",
             "Neutro2",
             "Mono/Macro1",
             "Mono/Macro2",
             "Lymphoid?"))

print(seu_zeb)
print(head(seu_zeb@meta.data))

# put the "cell types" into a slot with the same name as other datasets 
# cell_type, for downstream analysis compatibility
seu_zeb$cell_type <- seu_zeb$manually_annotated_ct

# add info on which column of the ensembl data frame to use based on Features 
seu_zeb@misc$ensembl_column_use <- "ENSDARG_ID" # zebrafish IDs

# test
stopifnot(!is.na(seu_zeb$cell_type))
stopifnot(!is.na(seu_zeb$manually_annotated_ct))

#-------------------------------------------------------------------------------

base::saveRDS(seu_zeb, snakemake@output[["zeb_seu"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# NAKED MOLE RAT

# this was readily available as seurat object, together with many other 
# mouse and NMR datasets. For now, stick to one main bone marrow dataset 
# (sorted)

#-------------------------------------------------------------------------------

Seurat_hgl_sorted_BM <- snakemake@input[["Seurat_hgl_sorted_BM"]]
Seurat_hgl_whole_BM <- snakemake@input[["Seurat_hgl_whole_BM"]]

# load objects into environment, then rename immediately
load(Seurat_hgl_sorted_BM)
seu_nmr_srt <- S.NOmp

load(Seurat_hgl_whole_BM)
seu_nmr_whl <- S.NOmp

#-------------------------------------------------------------------------------
# check some data
print(base::table(seu_nmr_srt$orig.ident))
print(base::table(seu_nmr_srt$tissue))
print(base::table(seu_nmr_srt$ident))
print(base::table(seu_nmr_srt$seurat_clusters))
print(base::table(seu_nmr_srt$celltype.combi))
print(base::table(seu_nmr_srt$celltype.combi.num))
print(base::table(seu_nmr_srt$sample))

#-------------------------------------------------------------------------------

# keep only cells derived from bone marrow, not blood

seu_nmr_srt <- seu_nmr_srt[,seu_nmr_srt$tissue == "marrow"]
                             
# remove cell types that are not in our dataset/conserved signatures
seu_nmr_srt <- seu_nmr_srt[,!seu_nmr_srt$celltype.combi %in% 
                             c("TC", "PC", "ERY", "MO", "DC", "PB-GC", "BC")]

# factorize cell types and put in correct slot for downstream analysis
seu_nmr_srt$celltype.combi <- factor(
  seu_nmr_srt$celltype.combi,
  levels = levels(seu_nmr_srt$celltype.combi)[
    levels(seu_nmr_srt$celltype.combi) %in% base::unique(
      seu_nmr_srt$celltype.combi)]
)

seu_nmr_srt$cell_type <- seu_nmr_srt$celltype.combi

print(base::table(seu_nmr_srt$cell_type))

stopifnot(!is.na(seu_nmr_srt$cell_type))

#-------------------------------------------------------------------------------

# add info on which column of the ensembl data frame to use based on Features 
seu_nmr_srt@misc$ensembl_column_use <- "HGLABER_SYMBOL" # NMR symbols

#-------------------------------------------------------------------------------

# remove umap coordinates and neighbors
seu_nmr_srt@neighbors <- list()
# keep original PCA coordinates but remove other reductions
seu_nmr_srt_pca <- seu_nmr_srt
seu_nmr_srt@reductions <- list()
seu_nmr_srt@reductions$pca_orig <- seu_nmr_srt_pca@reductions$pca

print("after removal")
print(seu_nmr_srt@reductions)

#-------------------------------------------------------------------------------

base::saveRDS(seu_nmr_srt, snakemake@output[["nmr_sorted_hspc"]])
base::saveRDS(seu_nmr_whl, snakemake@output[["nmr_whole_hspc"]])

utils::sessionInfo()
