#-------------------------------------------------------------------------------

library(DropletUtils)
library(scran)
library(bluster)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_11"]])
separate_fractions_clustering <- snakemake@params[["separate_fractions_clustering"]]

louvain <- clusterCells(sce, use.dimred="PCA", 
                        BLUSPARAM=NNGraphParam(cluster.fun="louvain"))

sce$louvain_clustering <- louvain

if(separate_fractions_clustering){
  print("separated fractions")
  
  # separate fractions for separate clustering 
  sce_hsc <- sce[,sce$Fraction_ID == "hsc"]
  sce_str <- sce[,sce$Fraction_ID == "str"]
  
  # cluster HSCs
  louvain_hsc <- clusterCells(sce_hsc, use.dimred="PCA", 
                              BLUSPARAM=NNGraphParam(cluster.fun="louvain"))
  sce_hsc$cluster_louvain <- as.numeric(unfactor(louvain_hsc))
  
  # cluster stromal cells
  louvain_str <- clusterCells(sce_str, use.dimred="PCA", 
                              BLUSPARAM=NNGraphParam(cluster.fun="louvain"))
  sce_str$cluster_louvain <- as.numeric(unfactor(louvain_str))  + 
    length(unique(sce_hsc$cluster_louvain))
  
  # put together again
  sce <- cbind(sce_hsc, sce_str)
  sce$cluster_louvain <- factor(sce$cluster_louvain, 
                               levels = sort(unique(sce$cluster_louvain)))
  sce$louvain_mode <- rep("fractions_separated", ncol(sce))
  
}else{
  
  louvain <- clusterCells(sce, use.dimred="PCA", 
                          BLUSPARAM=NNGraphParam(cluster.fun="louvain"))
  sce$cluster_louvain <- louvain
  sce$louvain_mode <- rep("fractions_together", ncol(sce))
  
}

print(sce)
saveRDS(sce, file = snakemake@output[["sce_12"]])