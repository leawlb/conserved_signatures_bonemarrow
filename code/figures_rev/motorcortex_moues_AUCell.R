library(Seurat)
library(AUCell)
library(tidyverse)
library(ggrastr)

# ----------------------------
# paths
# ----------------------------
proj_root <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data"

in_path <- file.path(proj_root, "metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS")

sce_dir <- file.path(proj_root, "scRNAseq/main_analysis/sce_objects/08_sce_brain")
nDEGs_path <- file.path(sce_dir, "01_list_nDEGs_all.rds")
cons_markers_path <- file.path(sce_dir, "02_marker_conserved_primates.rds")
out_rds <- file.path(sce_dir, "mouse_aucell_auc_meta.rds")

out_pdf <- file.path(proj_root, "manuscript1_rev/pdfs/motorcortex_mouse_AUCell.pdf")

# ----------------------------
# load + mouse subset
# ----------------------------
obj <- readRDS(in_path) |> UpdateSeuratObject()
mouse <- subset(obj, subset = orig.ident == "mouse" & subclass_label != "L6 IT Car3")

# ----------------------------
# make UMAP using existing PCA; no reclustering
# ----------------------------
mouse <- RunUMAP(mouse, reduction = "pca", dims = 1:50, n.components = 2)

# ----------------------------
# signatures: conserved = cons_markers ∩ nDEGs
# ----------------------------
nDEGs <- readRDS(nDEGs_path)
cons_markers <- readRDS(cons_markers_path)

conserved_signature <- lapply(names(cons_markers), function(ct) {
  intersect(cons_markers[[ct]], nDEGs[[ct]])
})
names(conserved_signature) <- names(cons_markers)
conserved_signature[["L6 IT Car3"]] <- NULL

# ----------------------------
# AUCell
# ----------------------------

expr <- LayerData(mouse[["RNA"]])
gene_sets_mouse <- lapply(conserved_signature, function(g) intersect(g, rownames(expr)))

rankings <- AUCell_buildRankings(expr, 
  plotStats = FALSE, 
  verbose = FALSE)
cellsAUC  <- AUCell_calcAUC(
  gene_sets_mouse,
  rankings,
  aucMaxRank = ceiling(0.1 * nrow(rankings))
)
auc_mat <- t(getAUC(cellsAUC))

mouse <- AddMetaData(mouse, auc_mat)
saveRDS(mouse, out_rds)

# ----------------------------
# plot
# ----------------------------

meta <- mouse@meta.data %>%
  dplyr::select(subclass_label, dplyr::all_of(colnames(auc_mat))) %>%
  tidyr::pivot_longer(
    cols = -subclass_label,
    names_to = "celltype",
    values_to = "score"
  )

order <- c("L2/3 IT", "L5 IT", "L6 IT", "L5 ET", "L5/6 NP", "L6 CT", "L6b")

violin <- ggplot(
  meta,
  aes(
    x = factor(subclass_label, levels = order),
    y = score
  )
) +
  facet_wrap(~ factor(celltype, levels = order), 
              ncol = 1, 
              scale = "free_y",
              strip.position = "right") +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none",
    strip.background = element_blank()
  ) +
  labs(
    x = NULL,
    y = "AUCell score"
  )


pdf(out_pdf, width = 7, height = 6)
  print(
    DimPlot(mouse, 
      group.by = "subclass_label", 
      reduction = "umap") & NoAxes()
  )
for (feat in colnames(auc_mat)) {
  p <- FeaturePlot(
    mouse,
    features = feat,
    reduction = "umap",
    raster = FALSE,
    order = TRUE
  ) +
    ggplot2::ggtitle(feat) &
    NoAxes()

  p$layers[[1]] <- ggrastr::rasterise(p$layers[[1]], dpi = 500)
  print(p)
}
print(violin)
dev.off()

