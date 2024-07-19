#-------------------------------------------------------------------------------
# Use seurats wilcoxon test-based method to find markers
# see https://www.biorxiv.org/content/10.1101/2022.05.09.490241v1

library(Seurat)
library(SeuratObject)
library(scran)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- base::readRDS(file = snakemake@input[["sce_input"]])
nr_hvgs <- snakemake@params[["nr_hvgs"]]

#-------------------------------------------------------------------------------
# convert to seurat
# use normalized logcounts
seurat <- Seurat::as.Seurat(sce, counts = "counts", data = "logcounts") 
SeuratObject::Idents(seurat) <- sce$cluster_louvain

#-------------------------------------------------------------------------------
# get hvgs for marker genes
hvgs <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(hvgs, n=nr_hvgs)

# find marker genes for each cluster
clustlist <- as.list(1:length(base::unique(sce$cluster_louvain)))

cluster_markers <- lapply(clustlist, function(x){
  
  markers <- Seurat::FindMarkers(seurat,
                                 slot = "data",
                                 test.use = "wilcox", 
                                 ident.1 = x, 
                                 features = hvgs)
  
  markers <- markers[base::order(base::abs(markers$avg_log2FC), 
                                 decreasing=TRUE),]
  
  markers$which_cluster <- base::rep(x, nrow(markers))
  
  return(markers)
})

#-------------------------------------------------------------------------------

base::saveRDS(cluster_markers, file = snakemake@output[["markers"]])

utils::sessionInfo()