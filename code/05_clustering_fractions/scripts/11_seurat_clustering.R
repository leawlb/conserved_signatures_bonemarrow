#-------------------------------------------------------------------------------

library(DropletUtils)
library(Seurat)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_09"]])
resolution <- snakemake@params[["resolution"]]

#-------------------------------------------------------------------------------
# Clustering

seurat <- as.Seurat(sce, counts = "counts", data = "logcounts") # logcounts (assay) only used if dims (reduction) = NULL

seurat <- FindNeighbors(seurat, dims = 1:10, reduction = "PCA") # batch corrected
seurat <- FindClusters(seurat, resolution = resolution)
sce$cluster_seurat <- Idents(seurat)
sce$cluster_seurat <- as.numeric(unfactor(sce$cluster_seurat))
sce$cluster_seurat <- sce$cluster_seurat + 1
sce$cluster_seurat <- factor(sce$cluster_seurat, 
                             levels = sort(unique(sce$cluster_seurat))) 
sce$cluster_seurat_resolution <- rep(resolution, ncol(sce))

print(sce)
saveRDS(sce, file = snakemake@output[["sce_11"]])