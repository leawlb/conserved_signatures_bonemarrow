
# get mouse stromal datasets

library(Seurat, quietly = TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# load dolgalev et al dataset used for cell type reference annotation
#seu_mus_str_ref <- base::readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/01_sce_prep/references_raw/bone-marrow-seurat.rds")
seu_mus_str_ref <- base::readRDS(snakemake@input[["seu_mus_str_ref"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# TIKHONOVA 

# get tikhonova  dataset only
seu_mus_tikhonova <- seu_mus_str_ref[,seu_mus_str_ref$orig.ident == "tikhonova"]

#-------------------------------------------------------------------------------
# check cell types

print(base::table(seu_mus_tikhonova$label))

# remove unknown cell type
ct_remove <- c("C")

seu_mus_tikhonova <- seu_mus_tikhonova[,which(!seu_mus_tikhonova$label == "C")]
print(base::table(seu_mus_tikhonova$label))

# factorize
seu_mus_tikhonova$label <- factor(
  seu_mus_tikhonova$label,
  levels = c("P1", "P2", "P3", "P4", "O1", "O2", "O3", "V1", "V2")
)

# put into cell_type slot
seu_mus_tikhonova$cell_type <- seu_mus_tikhonova$label 

print(base::table(seu_mus_tikhonova$cell_type))
stopifnot(!is.na(seu_mus_tikhonova$cell_type))

#-------------------------------------------------------------------------------

# check rownames
print(head(rownames(seu_mus_tikhonova)))

# add info on which ensembl column to use
seu_mus_tikhonova@misc$ensembl_column_use <- "MMUS_SYMBOL" # mouse gene symbols

#-------------------------------------------------------------------------------

# check assays

Seurat::DefaultAssay(seu_mus_tikhonova) <- "RNA"

print(seu_mus_tikhonova@assays)
print(seu_mus_tikhonova@assays$RNA)
print(seu_mus_tikhonova@assays$RNA@data[1:10,1:10])

# add logcounts into counts for consistency (with weinreb for example)
seu_mus_tikhonova@assays$RNA@counts <- seu_mus_tikhonova@assays$RNA@data

# remove "integrated" assay which has only 2000 features
seu_mus_tikhonova@assays <- list("RNA" = seu_mus_tikhonova@assays$RNA)
print(seu_mus_tikhonova@assays)

# add info on which slot to use
seu_mus_tikhonova@misc$data_use <- "logcounts" # only logcounts in assay

print("after removal")
print(seu_mus_tikhonova@assays$RNA)

#-------------------------------------------------------------------------------
# check reduced dims

print(seu_mus_tikhonova@reductions)

# remove umap coordinates and neighbors
seu_mus_tikhonova@neighbors <- list()
# keep original PCA coordinates but remove other reductions
seu_mus_tikhonova_pca <- seu_mus_tikhonova
seu_mus_tikhonova@reductions <- list()
#seu_mus_tikhonova@reductions$pca_orig <- seu_mus_tikhonova_pca@reductions$pca

print("after removal")
print(seu_mus_tikhonova@reductions)

#-------------------------------------------------------------------------------

base::saveRDS(seu_mus_tikhonova, snakemake@output[["mus_tik_stromal"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# BARYAWNO

# get baryawno dataset only
seu_mus_baryawno <- seu_mus_str_ref[,seu_mus_str_ref$orig.ident == "baryawno"]
print(dim(seu_mus_baryawno))
print(seu_mus_baryawno)

#-------------------------------------------------------------------------------
# check cell types

print(base::table(seu_mus_baryawno$label))

# factorize
seu_mus_baryawno$label <- factor(
  seu_mus_baryawno$label,
  levels = c(
    "Lepr-MSC-01",
    "OLC-1-07",
    "OLC-2-08",
    "Chondrocytes-progenitors-04",
    "Chondrocytes-13",
    "Chondrocytes-hypertrophic-02",
    "Chondrocytes-prehypertrophic-17",
    "Chondrocytes-proliferating-resting-10",
    "Fibroblasts-1-09",
    "Fibroblasts-2-15",
    "Fibroblasts-3-16",
    "Fibroblasts-4-03",
    "Fibroblasts-5-05",
    "EC-arterial-11",
    "EC-arteriolar-06",
    "EC-sinusoidal-00",
    "Pericytes-12"
  )
)

# put into cell_type slot
seu_mus_baryawno$cell_type <- seu_mus_baryawno$label 

print(base::table(seu_mus_baryawno$label))
stopifnot(!is.na(seu_mus_baryawno$label))

#-------------------------------------------------------------------------------

# check rownames
print(head(rownames(seu_mus_baryawno)))

# add info on which ensembl column to use
seu_mus_baryawno@misc$ensembl_column_use <- "MMUS_SYMBOL" # mouse gene symbols

#-------------------------------------------------------------------------------

# check assays

# no counts assay, no raw counts
Seurat::DefaultAssay(seu_mus_baryawno) <- "RNA"

print(seu_mus_baryawno@assays)
print(seu_mus_baryawno@assays$RNA)
print(seu_mus_baryawno@assays$RNA@data[1:10,1:10])

# add logcounts into counts for consistency (with weinreb for example)
seu_mus_baryawno@assays$RNA@counts <- seu_mus_baryawno@assays$RNA@data

# remove "integrated" assay which has only 2000 features
seu_mus_baryawno@assays <- list("RNA" = seu_mus_baryawno@assays$RNA)
print(seu_mus_baryawno@assays)

# add info on which slot to use
seu_mus_baryawno@misc$data_use <- "logcounts" # only logcounts in assay

print("after removal")
print(print(seu_mus_baryawno@assays$RNA))

#-------------------------------------------------------------------------------
# check reduced dims

print(seu_mus_baryawno@reductions)

# remove umap coordinates and neighbors
seu_mus_baryawno@neighbors <- list()
# keep original PCA coordinates but remove other reductions
seu_mus_baryawno_pca <- seu_mus_baryawno
seu_mus_baryawno@reductions <- list()
#seu_mus_baryawno@reductions$pca_orig <- seu_mus_baryawno_pca@reductions$pca

print(seu_mus_baryawno@reductions)

#-------------------------------------------------------------------------------

base::saveRDS(seu_mus_baryawno, snakemake@output[["mus_bar_stromal"]])
