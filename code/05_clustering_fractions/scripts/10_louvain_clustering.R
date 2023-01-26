#-------------------------------------------------------------------------------

library(DropletUtils)
library(scran)
library(bluster)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_09"]])
k_louvain <- snakemake@params[["k_louvain"]]

louvain <- clusterCells(sce, use.dimred="PCA", assay.type = NULL,
                        BLUSPARAM=NNGraphParam(cluster.fun="louvain", k = k_louvain))

sce$cluster_louvain <- louvain
sce$cluster_louvain <- factor(sce$cluster_louvain, 
                              levels = sort(unique(sce$cluster_louvain))) 
sce$k_louvain <- rep(k_louvain, ncol(sce))

print(sce)
saveRDS(sce, file = snakemake@output[["sce_10"]])