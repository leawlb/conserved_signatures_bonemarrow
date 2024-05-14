#-------------------------------------------------------------------------------
# Use seurats wilcoxon test-based method to find markers in each cell type
# see https://www.biorxiv.org/content/10.1101/2022.05.09.490241v1

library(Seurat)
library(SeuratObject)
library(scran)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# convert to seurat
seurat <- Seurat::as.Seurat(sce, counts = "counts", data = "logcounts",) 
SeuratObject::Idents(seurat) <- sce$celltypes
print(unique(sce$celltypes))

#-------------------------------------------------------------------------------
# get hvgs for marker genes
hvgs <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(hvgs, n=nr_hvgs)

# find marker genes for each celltype
ctlist <- as.list(unique(sce$celltypes))
ct_markers <- lapply(ctlist, function(x){
  markers <- Seurat::FindMarkers(seurat,
                                 test.use = "wilcox", ident.1 = x, 
                                 features = hvgs)
  markers <- markers[order(abs(markers$avg_log2FC), decreasing=TRUE),]
  markers$which_cluster <- rep(x, nrow(markers))
  return(markers)
})

#-------------------------------------------------------------------------------
saveRDS(ct_markers, file = snakemake@output[["markers"]])

sessionInfo()