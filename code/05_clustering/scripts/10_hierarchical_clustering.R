#-------------------------------------------------------------------------------

library(DropletUtils)
library(bluster)
source(file = snakemake@params[["sce_functions"]])

#-------------------------------------------------------------------------------


sce <- readRDS(file = snakemake@input[["sce_09"]])
remove_mature_cells_clustering <- snakemake@params[["remove_mature_cells_clustering"]]
separate_fractions_clustering <- snakemake@params[["separate_fractions_clustering"]]
number_k <- snakemake@params[["number_k"]]

# remove mature cells carefully (TCs and BCs only) based on ref annotations
sce <- get_mature_cts_all(sce)
if(remove_mature_cells_clustering){
  sce <- sce[,sce$mature_cells == "FALSE"]
}
# only a small number of cells is affected by removal of TCs and BCs,
# but cast and caroli contain large numbers of putative contamination
  
# separate fractions for separate clustering 
sce_hsc <- sce[,sce$Fraction_ID == "hsc"]
sce_str <- sce[,sce$Fraction_ID == "str"]

# find contamination but do not remove it (save it for later)
sce_hsc <- find_contamination_ref_all(sce_hsc, input = "hsc")
sce_str <- find_contamination_ref_all(sce_str, input = "stromal")

#-------------------------------------------------------------------------------

# extract correct number of ks used for cutting the dendrogram
species_curr <- sce$Species_ID[1]

number_k_hsc <- number_k[[species_curr]]$hsc
number_k_str <- number_k[[species_curr]]$str
number_k_both <- number_k[[species_curr]]$both
print(number_k_hsc)
print(number_k_str)

#-------------------------------------------------------------------------------

if(separate_fractions_clustering){
  print("separated fractions")
  
  # HSC clustering
  mat_hsc <- reducedDim(sce_hsc, "PCA") # PCA is batch corrected
  mat_hsc <- scale(mat_hsc)
  dist_hsc <- dist(mat_hsc, method = 'euclidean')
  ward_hsc <- hclust(dist_hsc, method = "ward.D2")

  cut_hsc <- cutree(ward_hsc, k = number_k_hsc)
  sce_hsc$cluster_hierarchical <- cut_hsc

  # Stromal clustering
  mat_str <- reducedDim(sce_str, "PCA")
  mat_str <- scale(mat_str)
  dist_str <- dist(mat_str, method = 'euclidean')
  ward_str <- hclust(dist_str, method = "ward.D2")

  cut_str <- cutree(ward_str, k = number_k_str)
  sce_str$cluster_hierarchical <- cut_str

  # add together
  sce_str$cluster_hierarchical <- sce_str$cluster_hierarchical + number_k_hsc
  sce <- cbind(sce_hsc, sce_str)

  sce$cluster_hierarchical <- factor(sce$cluster_hierarchical, 
                                     levels = sort(unique(sce$cluster_hierarchical)))
  
  sce$hierarchical_mode <- rep("fractions_separated", ncol(sce))

}else{
  
  mat <- reducedDim(sce, "PCA") # PCA is batch corrected
  mat <- scale(mat)
  dist <- dist(mat, method = 'euclidean')
  ward <- hclust(dist, method = "ward.D2")
  
  cut <- cutree(ward, k = number_k_both)
  
  sce$cluster_hierarchical <- cut
  sce$cluster_hierarchical <- factor(sce$cluster_hierarchical, 
                                     levels = sort(unique(sce$cluster_hierarchical)))
  sce$hierarchical_mode <- rep("fractions_together", ncol(sce))
  
}

print(sce)
saveRDS(sce, file = snakemake@output[["sce_10"]])
