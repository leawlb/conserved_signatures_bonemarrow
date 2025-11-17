library(Seurat)
library(tidyverse)
library(ggrastr)

mouse_lemur_harmony <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/mouse_lemur_Ezran2025/02_mouse_lemur_harmony.rds")

celltypes <- c("Adipo.CAR", 
                "Art.Cap.EC",
                "BM.prog",
                "Cycling",
                "Early.MPP",
                "Early.Myeloid",
                "Erythroid",
                "Fibro.Chondro",
                "Fibroblast",
                "HSC",
                "Lymphoid",
                "Mk.prog",
                "Mk.Ery.prog",
                "Mono.prog",
                "Neutro.prog",
                "Osteo",
                "Peri.SMC",
                "Sinusoid.EC",
                "Sinusoid.Ven.EC")

pdf(
  file   = "code/figures/revision/mouse_leumur_AUCell_weights.pdf",
  height = 7, width = 7
)
for (feat in celltypes) {
  p <- FeaturePlot(
    mouse_lemur_harmony,
    features  = feat,
    reduction = "umap_harmony",
    raster    = FALSE,   # let FeaturePlot create normal vector points
    order     = TRUE
  ) + ggplot2::ggtitle(feat) & NoAxes()

  ## Rasterize the first layer (the points) at high DPI
  p$layers[[1]] <- ggrastr::rasterise(p$layers[[1]], dpi = 500)

  print(p)
}
dev.off()


DimPlot(mouse_lemur_harmony,
  reduction = "umap_harmony",
  label     = TRUE,
  repel     = TRUE,
  raster    = FALSE
) & NoAxes() & NoLegend()

## 21 = HSCs
## 26 & 29 = stromal/EC


