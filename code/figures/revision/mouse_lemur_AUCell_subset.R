library(Seurat)
library(harmony)

mouse_lemur_harmony <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/mouse_lemur_Ezran2025/02_mouse_lemur_harmony.rds")

HSC <- subset(mouse_lemur_harmony,
              idents = 21)
DefaultAssay(HSC) <- "RNA" 

HSC <- HSC |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000) |>
  ScaleData() |>
  RunPCA(npcs = 50, verbose = FALSE)

HSC <- harmony::RunHarmony(
  object         = HSC,
  group.by.vars  = "tissue", 
  reduction.use  = "pca",
  dims.use       = 1:50,
  reduction.save = "harmony",
  project.dim    = FALSE,
  verbose        = TRUE
)

HSC <- HSC |>
  FindNeighbors(reduction = "harmony", dims = 1:50) |>
  RunUMAP(
    reduction      = "harmony",
    dims           = 1:50,
    reduction.name = "umap_harmony",
    n.neighbors    = 200,
    min.dist       = 1
  ) |>
  FindClusters(resolution = 0.5)


pdf(
  file   = "code/figures/revision/mouse_leumur_AUCell_HSC.pdf",
  height = 3, width = 3)
FeaturePlot(
    HSC,
    features  = "HSC",
    reduction = "umap_harmony",
    raster    = FALSE,
    order     = TRUE
) & NoAxes()
dev.off()





stromal <- subset(mouse_lemur_harmony,
                  idents = c(26, 29))
DefaultAssay(stromal) <- "RNA" 

stromal <- stromal |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000) |>
  ScaleData() |>
  RunPCA(npcs = 50, verbose = FALSE)

stromal <- harmony::RunHarmony(
  object         = stromal,
  group.by.vars  = "tissue",
  reduction.use  = "pca",
  dims.use       = 1:50,
  reduction.save = "harmony",
  project.dim    = FALSE,
  verbose        = TRUE
)

stromal <- stromal |>
  FindNeighbors(reduction = "harmony", dims = 1:50) |>
  RunUMAP(
    reduction      = "harmony",
    dims           = 1:10,
    reduction.name = "umap_harmony",
    n.neighbors    = 200,
    min.dist       = 1
  ) |>
  FindClusters(resolution = 0.5)


celltypes <- c("Adipo.CAR", 
                "Art.Cap.EC",
                "Cycling",
                "Fibro.Chondro",
                "Fibroblast",
                "Osteo",
                "Peri.SMC",
                "Sinusoid.EC",
                "Sinusoid.Ven.EC")

pdf(file = "code/figures/revision/mouse_leumur_AUCell_stromal.pdf",
  height = 4, width = 4)
for (feat in celltypes) {
  p <- FeaturePlot(
    stromal,
    features  = feat,
    reduction = "umap_harmony",
    raster    = FALSE,
    order     = TRUE
  ) + ggplot2::ggtitle(feat) & NoAxes()

  ## Rasterize the first layer (the points) at high DPI
  p$layers[[1]] <- ggrastr::rasterise(p$layers[[1]], dpi = 500)

  print(p)
}
dev.off()
