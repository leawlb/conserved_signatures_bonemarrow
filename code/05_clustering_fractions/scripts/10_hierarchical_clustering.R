#-------------------------------------------------------------------------------

library(DropletUtils)
library(bluster)
source(file = snakemake@params[["sce_functions"]])
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_09"]])
number_k <- snakemake@params[["number_k"]]

# extract correct number of ks used for cutting the dendrogram
fraction_curr <- sce$Fraction_ID[1]
number_k_curr <- number_k[[fraction_curr]]
print(number_k_curr)

#-------------------------------------------------------------------------------

mat <- reducedDim(sce, "PCA") # PCA is batch corrected
mat <- scale(mat)
dist <- dist(mat, method = 'euclidean')
ward <- hclust(dist, method = "ward.D2")
  
cut <- cutree(ward, k = number_k_both)
  
sce$cluster_hierarchical <- cut
sce$cluster_hierarchical <- factor(sce$cluster_hierarchical, 
                                   levels = sort(unique(sce$cluster_hierarchical)))

print(sce)
saveRDS(sce, file = snakemake@output[["sce_10"]])
