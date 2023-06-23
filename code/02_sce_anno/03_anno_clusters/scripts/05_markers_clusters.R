#-------------------------------------------------------------------------------

library(DropletUtils)
library(Seurat)
library(scran)
set.seed(37)

# Use seurats wilcoxon test-based method to find markers
# see https://www.biorxiv.org/content/10.1101/2022.05.09.490241v1

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# convert to seurat
seurat <- as.Seurat(sce, counts = "counts", data = "logcounts",) 
Idents(seurat) <- sce$cluster_louvain

# get hvgs for marker genes
hvgs <- modelGeneVar(sce)
hvgs <- getTopHVGs(hvgs, n=nr_hvgs)

# find marker genes for each cluster
clustlist <- as.list(1:length(unique(sce$cluster_louvain)))
cluster_markers <- lapply(clustlist, function(x){
  markers <- FindMarkers(seurat, test.use = "wilcox", ident.1 = x, 
                         features = hvgs)
  markers <- markers[order(abs(markers$avg_log2FC), decreasing=TRUE),]
  markers$which_cluster <- rep(x, nrow(markers))
  return(markers)
})

saveRDS(cluster_markers, file = snakemake@output[["markers"]])
