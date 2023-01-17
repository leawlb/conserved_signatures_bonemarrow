#-------------------------------------------------------------------------------

library(DropletUtils)
library(bluster)
source(file = snakemake@params[["sce_functions"]])
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_09"]])
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
  sce <- sce[,sample_numbers]
  
  # recalculate PCA from batch corrected values
  hvgs <- modelGeneVar(sce)
  hvgs <- getTopHVGs(hvgs, n=nr_hvgs)
  if(sce$Correction_method == "FastMNN"){
    sce <- runPCA(sce, ncomponents=25, subset_row = hvgs, 
                  exprs_values = "reconstructed")
  }else if(sce$Correction_method == "Seurat"){
    sce <- runPCA(sce, ncomponents=25, subset_row = hvgs, 
                  exprs_values = "corrected") 
  }
}

#-------------------------------------------------------------------------------

mat <- reducedDim(sce, "PCA") # PCA is batch corrected
mat <- scale(mat[,1:10]) 
dist <- dist(mat, method = 'euclidean')
ward <- hclust(dist, method = "ward.D2")
  
cut <- cutree(ward, k = number_k_curr)
  
sce$cluster_hierarchical <- cut
sce$cluster_hierarchical <- factor(sce$cluster_hierarchical, 
                                   levels = sort(unique(sce$cluster_hierarchical)))

print(sce)
saveRDS(sce, file = snakemake@output[["sce_10"]])
