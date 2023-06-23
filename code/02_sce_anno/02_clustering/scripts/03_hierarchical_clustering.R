#-------------------------------------------------------------------------------

library(DropletUtils)
library(bluster)

source("../../source/sce_functions.R")
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
number_k <- snakemake@params[["number_k"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]
size_subs_hscs <- snakemake@params[["size_subs_hscs"]]

# extract correct number of ks used for cutting the dendrogram from config
fraction_curr <- sce$Fraction_ID[1]
number_k_curr <- number_k[[fraction_curr]]
print(number_k_curr)

# sub sample HSC object since it's too big
if(fraction_curr == "hsc"){
  RNGkind("L'Ecuyer-CMRG") # for random number generation
  set.seed(37)
  sample_numbers <- sample(x = 1:ncol(sce), size = size_subs_hscs)
  print(sample_numbers[123])
  sce <- sce[,sample_numbers]
}

#-------------------------------------------------------------------------------

set.seed(37)
mat <- reducedDim(sce, "PCA") # PCA is batch corrected
mat <- scale(mat[,1:10]) 
dist <- dist(mat, method = 'euclidean')
ward <- hclust(dist, method = "ward.D2")
  
cut <- cutree(ward, k = number_k_curr)
  
sce$cluster_hierarchical <- cut
sce$cluster_hierarchical <- factor(
  sce$cluster_hierarchical, levels = sort(unique(sce$cluster_hierarchical)))

print(sce)
saveRDS(sce, file = snakemake@output[["sce_03"]])