library(Seurat)
library(harmony)
library(future)

options(future.globals.maxSize = 30 * 1024^3)

setwd("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/mouse_lemur_Ezran2025")

mouse_lemur_rawcounts <- readRDS("exports/combined_seurat_obj.rds")

# DimPlot(
#   mouse_lemur_rawcounts,
#   reduction = "umap",
#   group.by  = "cell_ontology_class_v1",
#   split.by  = "tissue",
#   label     = TRUE,
#   raster    = TRUE,
#   ncol = 1
# ) & NoAxes() & NoLegend()


mouse_lemur_harmony <- mouse_lemur_rawcounts |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000) |>
  ScaleData() |>
  RunPCA(npcs = 50, verbose = FALSE)

mouse_lemur_harmony <- harmony::RunHarmony(
  object        = mouse_lemur_harmony,
  group.by.vars = c("tissue"),
  reduction.use = "pca",
  dims.use      = 1:50,
  reduction.save= "harmony",
  project.dim   = FALSE,
  verbose       = TRUE
)

# Graph + UMAP on Harmony embeddings
mouse_lemur_harmony <- mouse_lemur_harmony |>
  FindNeighbors(reduction = "harmony", 
    dims = 1:50) |>
  RunUMAP(reduction = "harmony", 
    dims = 1:50, 
    reduction.name = "umap_harmony") |>
  FindClusters(resolution = 0.5)

DimPlot(mouse_lemur_harmony, 
  reduction = "umap_harmony", 
  group.by = "tissue", 
  raster = TRUE) & NoAxes()

DimPlot(mouse_lemur_harmony,
  reduction = "umap_harmony",
  group.by  = "free_annotation_v1",
  label     = TRUE,
  repel     = TRUE,
  raster    = FALSE
) & NoAxes() & NoLegend()

DimPlot(mouse_lemur_harmony,
  reduction = "umap_harmony",
  group.by  = "compartment_v1",
  label     = FALSE,
  raster    = FALSE
) & NoAxes()

saveRDS(mouse_lemur_harmony, 
  file = "01_mouse_lemur_harmony.rds")

