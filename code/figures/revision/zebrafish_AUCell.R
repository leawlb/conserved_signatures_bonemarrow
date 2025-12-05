library(Seurat)
library(tidyverse)
library(ggrastr)

zebrafish_HSC <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/zebrafish_Athanasiadis2017/03_zebrafish_new.rds")

celltypes <- c( "Cycling",
                "HSC",
                "Early.MPP",
                "Mk.prog",
                "BM.prog",
                "Mk.Ery.prog",
                "Erythroid",
                "Early.Myeloid",
                "Mono.prog",
                "Neutro.prog",
                "Lymphoid")

meta <- zebrafish_HSC@meta.data[,c("Component1", "Component2", "Type", celltypes)]

pdf(
  file   = "code/figures/revision/zebrafish_AUCell_weights_new.pdf",
  height = 2, width = 3
)

ggplot(meta,
aes(x = Component1,
    y = Component2,
    color = Type)) +
      geom_point(size = 0.5) +
      theme_void()

for (feat in celltypes) {
  meta_sorted <- meta |>
    arrange(.data[[feat]])   # sort by the feature column, low → high

  p <- ggplot(
      meta_sorted,
      aes(
        x     = Component1,
        y     = Component2,
        color = .data[[feat]]
      )
    ) +
    geom_point(size = 0.5) +
    theme_void() +
    scale_color_gradient(
    low    = "grey90",
    high   = "blue",
    limits = c(0, 0.8),    # lower bound fixed at 0; upper bound free
    oob    = scales::squish) +
    labs(color = NULL) +
    ggtitle(feat)

  print(p)

}
dev.off()

plot_meta <- gather(meta, celltype, score, -Component1, -Component2, -Type)


pdf(
  file   = "code/figures/revision/zebrafish_AUCell_weights_new_clusters.pdf",
  height = 3, width = 5)
ggplot(filter(plot_meta, celltype %in% c("Erythroid", "HSC", "Neutro.prog", "Mono.prog", "Mk.prog")),
      aes(x = factor(Type, levels = c("HSPC", "Neutrophils", "Monocytes", "Erythrocytes", "Thrombocytes")),
        y = score,
        color = Type)) +
  facet_wrap(~factor(celltype, levels = c("HSC", "Neutro.prog", "Mono.prog", "Erythroid", "Mk.prog")), nrow = 1) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, h = 1, v = 0.5),
        legend.position = "none") +
  labs(x = NULL, y = "AUCell score")
dev.off()
