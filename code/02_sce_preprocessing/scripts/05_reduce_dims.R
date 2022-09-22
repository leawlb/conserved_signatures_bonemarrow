#-------------------------------------------------------------------------------
# author: Amy Danson, Lea Wölbert
# date: 2022-09-21
# extract highly variable genes and reduce dimensions for QC purposes

library(scater, quietly = TRUE) 
library(scran, quietly = TRUE) 

sce <- readRDS(file = snakemake@input[["sce_04"]])

# extract hvgs for subsequent dimensionality reduction at sample level
genevar <- modelGeneVar(sce)
hvg <- getTopHVGs(genevar, n=2000)

# reduce dimensions
sce <- runPCA(sce, ncomponents=25, subset_row = hvg) 
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE,  subset_row = hvg)

# save as new object
saveRDS(sce, file = snakemake@output[["sce_05"]])
