#-------------------------------------------------------------------------------

library(DropletUtils)
library(bluster)
source(file = snakemake@params[["sce_functions"]])

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_07"]])
number_k <- snakemake@params[["number_k"]]
print(number_k)

# remove mature cells carfully (TCs and BCs only) based on ref annotations
sce <- get_mature_cts_all(sce)
sce <- sce[,sce$mature_cells == "FALSE"]

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
species_curr <- sce$Species[1]
number_k = number_k[species_curr %in% number_k]

number_k_hsc = number_k["hsc" %in% number_k]
number_k_str = number_k["str" %in% number_k]

#-------------------------------------------------------------------------------

# HSC clustering
mat_hsc <- reducedDim(sce_hsc, "PCA")
mat_hsc <- scale(mat_hsc)
dist_hsc <- dist(mat_hsc, method = 'euclidean')
ward_hsc <- hclust(dist_hsc, method = "ward.D2")

cut_hsc <- cutree(ward_hsc, k = number_k_hsc)
sce_hsc$hclust <- cut_hsc

# Stromal clustering
mat_str <- reducedDim(sce_str, "PCA")
mat_str <- scale(mat_str)
dist_str <- dist(mat_str, method = 'euclidean')
ward_str <- hclust(dist_str, method = "ward.D2")

cut_str <- cutree(ward_str, k = number_k_str)
sce_str$hclust <- cut_str

sce_str$hclust <- sce_str$hclust + number_k_hsc

sce <- cbind(sce_hsc, sce_str)


sce$hclust <- factor(sce$hclust, levels = unique(sce$hclust))



# could also try HGC and scran clusterCells, Seurat