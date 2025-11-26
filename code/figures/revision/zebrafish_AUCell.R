library(Seurat)
library(tidyverse)
library(ggrastr)

zebrafish_HSC <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/zebrafish_Athanasiadis2017/03_zebrafish.rds")

celltypes <- c("BM.prog", 
                "Cycling",
                "Early.MPP",
                "Early.Myeloid",
                "Erythroid",
                "HSC",
                "Lymphoid",
                "Mk.prog",
                "Mk.Ery.prog",
                "Mono.prog",
                "Neutro.prog")

pdf(
  file   = "code/figures/revision/zebrafish_AUCell_weights.pdf",
  height = 4, width = 4
)
for (feat in celltypes) {
  p <- FeaturePlot(
    zebrafish_HSC,
    features  = feat,
    reduction = "umap_pca",
    raster    = FALSE,
    order     = TRUE) + 
    ggplot2::ggtitle(feat) +
    ggplot2::scale_color_gradient(limits = c(0, 0.6),
    oob = scales::squish,
    low = "grey80", high = "blue") & NoAxes()

  ## Rasterize the first layer (the points) at high DPI
  p$layers[[1]] <- ggrastr::rasterise(p$layers[[1]], dpi = 500)

  print(p)
}
dev.off()


DimPlot(zebrafish_HSC,
  reduction = "umap_pca",
  group.by = "cell_type", 
  label     = TRUE,
  repel     = TRUE,
  raster    = FALSE
) & NoAxes() & NoLegend()



