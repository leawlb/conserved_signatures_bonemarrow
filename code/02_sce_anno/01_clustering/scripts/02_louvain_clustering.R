#-------------------------------------------------------------------------------
# louvain clustering based on corrected PCs

library(scran, quietly = TRUE)
library(igraph, quietly = TRUE)
#library(bluster, quietly = TRUE) # for NNGraphParam()
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
k_graph_list <- snakemake@params[["k_graph_list"]]
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]
resolution_louvain <- resolution_louvain_list[[fraction_curr]]
print(k_graph)
print(resolution_louvain)

graph <- scran::buildSNNGraph(sce, 
                              k = k_graph, 
                              use.dimred = "PCA", # corrected
                              type = "rank")

clust <- igraph::cluster_louvain(graph, 
                                 resolution = resolution_louvain)

sce$cluster_louvain <- as.character(clust$membership)
sce$cluster_louvain <- factor(sce$cluster_louvain, 
                              levels = sort(unique(sce$cluster_louvain))) 
sce$k_louvain <- rep(k_graph, ncol(sce))
sce$resolution_louvain <- rep(resolution_louvain, ncol(sce))

print(sce)
saveRDS(sce, file = snakemake@output[["sce_output"]])

sessionInfo()