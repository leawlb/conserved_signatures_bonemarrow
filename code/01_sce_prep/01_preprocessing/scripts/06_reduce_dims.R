#-------------------------------------------------------------------------------
# authors: Amy Danson, Lea Wölbert
# extract highly variable genes and reduce dimensions for QC purposes

library(scater, quietly = TRUE) 
library(scran, quietly = TRUE) 
set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_input"]])
nr_hvgs <- snakemake@params[["nr_hvgs"]]

genevar <- modelGeneVar(sce)
hvg <- getTopHVGs(genevar, n=nr_hvgs)

# reduce dimensions
set.seed(37)
sce <- runPCA(sce, ncomponents=25, subset_row = hvg) 
set.seed(37)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE, subset_row = hvg)

saveRDS(sce, file = snakemake@output[["sce_output"]])