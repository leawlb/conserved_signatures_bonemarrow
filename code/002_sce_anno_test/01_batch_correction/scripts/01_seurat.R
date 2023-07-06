#-------------------------------------------------------------------------------

# instructions from https://satijalab.org/seurat/archive/v3.0/integration.html

library(DropletUtils, quietly = TRUE)
library(Seurat, quietly = TRUE)
source(file = "../../source/sce_functions.R")
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]]) #doesn't need renormalization
batch_use <- snakemake@params[["batch_use"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]
nr_hvgs_batch_correction <- snakemake@params[["nr_hvgs_batch_correction"]]
seurat_reduction <- snakemake@params[["seurat_reduction"]]

#-------------------------------------------------------------------------------
# proceed with seurat functionalities
seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")

seurat_list <- SplitObject(seurat, split.by = batch_use)

for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = nr_hvgs_batch_correction, 
                                           verbose = FALSE)
}

if(seurat_reduction == "rpca"){
  features <- SelectIntegrationFeatures(object.list = seurat_list)
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
}

seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30,
                                         anchor.features = nr_hvgs_batch_correction,
                                         reduction = seurat_reduction)
print("starting IntegrateData")
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)

sce_seu <- as.SingleCellExperiment(seurat_integrated)

#-------------------------------------------------------------------------------

# rename the already existent Reduced Dimensions
reducedDimNames(sce)[grep("PCA", reducedDimNames(sce))] <- "PCA_before_BC"
reducedDimNames(sce)[grep("UMAP", reducedDimNames(sce))] <- "UMAP_before_BC"

#-------------------------------------------------------------------------------

# transfer to sce

# subset SCE to remaining genes and add corrected values
sce <- sce[rownames(sce) %in% rownames(sce_seu),]
assays(sce)$corrected <- assays(sce_seu)$logcounts[
  match(rownames(assays(sce)$logcounts), rownames(assays(sce_seu)$logcounts)),
  match(colnames(assays(sce)$logcounts), colnames(assays(sce_seu)$logcounts))]

# calculate new PC and UMAP coordinates based on corrected values
hvgs <- modelGeneVar(sce, assay.type = "corrected")
hvgs <- getTopHVGs(hvgs, n=nr_hvgs)
sce <- runPCA(sce, ncomponents=25, subset_row = hvgs, 
              exprs_values = "corrected") 
sce <- runUMAP(sce, dimred = 'PCA', exprs_values = "corrected",
               external_neighbors=TRUE, subset_row = hvgs)

colData(sce)$Correction_method <- "Seurat"

print(sce)
saveRDS(sce, file = snakemake@output[["sce_output"]]) 