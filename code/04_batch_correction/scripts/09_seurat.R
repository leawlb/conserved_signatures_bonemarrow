#-------------------------------------------------------------------------------

# instructions from https://satijalab.org/seurat/archive/v3.0/integration.html

library(DropletUtils, quietly = TRUE)
library(Seurat)

sce <- readRDS(file = snakemake@input[["sce_07"]]) # doesn't need renormalization
batch_use <- snakemake@params[["batch_use"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]
  
seurat <- as.Seurat(sce[unique(rownames(sce)),],
                    counts = "counts", data = "logcounts")

seurat_list <- SplitObject(seurat, split.by = batch_use)

for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = nr_hvgs, verbose = FALSE)
}

seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)

sce <- as.SingleCellExperiment(seurat_integrated)

saveRDS(sce, file = snakemake@outpur[["sce_09"]]) 
# this has 40 removed rows because of duplication compared to other means of BC


