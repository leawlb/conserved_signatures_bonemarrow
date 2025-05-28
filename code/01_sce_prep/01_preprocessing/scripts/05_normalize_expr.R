#-------------------------------------------------------------------------------
# authors: Amy Danson, Lea Wölbert
# normalize and log-transform expression on the level of individual samples 

RNGkind("L'Ecuyer-CMRG") 

library(scran, quietly = TRUE) 
set.seed(37)

sce <- base::readRDS(file = snakemake@input[["sce_input"]])

# quick clustering to accelerate calculation
quick_clust <- scran::quickCluster(sce) 

# add sum factors (size factor equivalent) to SCE, normalize, log-transform
sce <- scran::computeSumFactors(sce, cluster = quick_clust)
sce <- scuttle::logNormCounts(sce) 

base::saveRDS(sce, file = snakemake@output[["sce_output"]])

utils::sessionInfo()