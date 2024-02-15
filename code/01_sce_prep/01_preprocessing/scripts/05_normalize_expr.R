#-------------------------------------------------------------------------------
# authors: Amy Danson, Lea Wölbert
# normalize and log-transform expression on the level of individual samples 

library(SingleCellExperiment, quietly = TRUE) 
library(scran, quietly = TRUE) 
set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_input"]])

quick_clust <- quickCluster(sce) # quick clustering to accelerate calculation

# must be COMPUTE not CALCULATE SumFactors to add directly to SCE
sce <- computeSumFactors(sce, cluster = quick_clust)
sce <- logNormCounts(sce) # log-transformation

saveRDS(sce, file = snakemake@output[["sce_output"]])