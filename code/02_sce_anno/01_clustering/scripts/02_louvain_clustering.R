#-------------------------------------------------------------------------------
# louvain clustering based on corrected PCs

library(scran, quietly = TRUE)
library(igraph, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# objects and params 
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

k_graph_list <- snakemake@params[["k_graph_list"]]
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]
resolution_louvain <- resolution_louvain_list[[fraction_curr]]

print(k_graph)
print(resolution_louvain)

#-------------------------------------------------------------------------------
# get a graph of nearest neighbors
graph <- scran::buildSNNGraph(sce, 
                              k = k_graph, 
                              use.dimred = "PCA", # corrected
                              type = "rank")

# community detection
clust <- igraph::cluster_louvain(graph, 
                                 resolution = resolution_louvain)

#-------------------------------------------------------------------------------
# add to sce object 
sce$cluster_louvain <- clust$membership

sce$cluster_louvain <- base::factor(
  sce$cluster_louvain, 
  levels = base::sort(base::unique(sce$cluster_louvain)))


sce$k_graph <- base::rep(k_graph, ncol(sce))
sce$resolution_louvain <- base::rep(resolution_louvain, ncol(sce))

print(sce)

#-------------------------------------------------------------------------------
base::saveRDS(sce, file = snakemake@output[["sce_output"]])

utils::sessionInfo()