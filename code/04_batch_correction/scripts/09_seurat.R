#-------------------------------------------------------------------------------

# instructions from https://satijalab.org/seurat/archive/v3.0/integration.html

library(DropletUtils, quietly = TRUE)
library(Seurat)

sce <- readRDS(file = snakemake@input[["sce_07"]]) # doesn't need renormalization
batch_use <- snakemake@params[["batch_use"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]
nr_hvgs_batch_correction <- snakemake@params[["nr_hvgs_batch_correction"]]

print(sce)
print(batch_use)
print(nr_hvgs_batch_correction)

source(file = snakemake@params[["sce_functions"]])

# rename the already existent Reduced Dimensions
reducedDimNames(sce)[grep("PCA", reducedDimNames(sce))] <- "PCA_before"
reducedDimNames(sce)[grep("UMAP", reducedDimNames(sce))] <- "UMAP_before"

seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")

seurat_list <- SplitObject(seurat, split.by = batch_use)

for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = nr_hvgs_batch_correction, 
                                           verbose = FALSE)
}

seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                         anchor.features = nr_hvgs_batch_correction)
print("starting IntegrateData")
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)

#sce <- as.SingleCellExperiment(seurat_integrated)
sce_seu <- as.SingleCellExperiment(seurat_integrated)

# subset SCE to remaining genes and add corrected values
sce <- sce[rownames(sce) %in% rownames(sce_seu),]
assays(sce, withDimnames=FALSE)$corrected <- assays(sce_seu)$logcounts
# calculate new PC and UMAP coordinates based on corrected values
#sce <- reduce_dims(sce, nr_hvgs = nr_hvgs) # new values stored in logcounts
hvgs <- modelGeneVar(sce, assay.type = "corrected")
hvgs <- getTopHVGs(hvgs, n=nr_hvgs)
sce <- runPCA(sce, ncomponents=25, subset_row = hvgs, 
              exprs_values = "corrected") 
sce <- runUMAP(sce, dimred = 'PCA', exprs_values = "corrected",
               external_neighbors=TRUE, subset_row = hvgs)
saveRDS(sce, file = snakemake@output[["sce_09"]]) 

