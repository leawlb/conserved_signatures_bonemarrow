#-------------------------------------------------------------------------------

library(DropletUtils)
library(scran)
library(bluster)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_11"]])

louvain <- clusterCells(sce, use.dimred="PCA", assay.type = NULL,
                        BLUSPARAM=NNGraphParam(cluster.fun="louvain"))

sce$louvain_clustering <- louvain

print(sce)
saveRDS(sce, file = snakemake@output[["sce_12"]])