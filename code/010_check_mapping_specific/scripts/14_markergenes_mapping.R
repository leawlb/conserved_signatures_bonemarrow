#-------------------------------------------------------------------------------

library(biomaRt)
library(Seurat) 
library(scran) 
set.seed(37)

#-------------------------------------------------------------------------------

sce_tog <- readRDS(file = snakemake@input[["sce_tog"]])

nr_hvgs <- snakemake@params[["nr_hvgs"]]

# convert to seurat
seurat <- as.Seurat(sce_tog, counts = "counts", data = "logcounts",) 
Idents(seurat) <- sce_tog$Genome

# get hvgs for marker genes
hvgs <- modelGeneVar(sce_tog)
hvgs <- getTopHVGs(hvgs, n=nr_hvgs)

# find marker genes for each cluster
genomes <- unique(sce_tog$Genome)
clustlist <- as.list(genomes)
cluster_markers <- lapply(clustlist, function(x){
  markers <- FindMarkers(seurat, test.use = "wilcox", ident.1 = x, 
                         features = hvgs)
  markers <- markers[order(abs(markers$avg_log2FC), decreasing=TRUE),]
  markers$which_cluster <- rep(x, nrow(markers))
  return(markers)
})

names(cluster_markers) <- genomes
markergenes_og <- cluster_markers[["OneGenome"]]

saveRDS(markergenes_og, file = snakemake@output[["markergenes_mapping"]])