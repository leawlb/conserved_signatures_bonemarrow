#-------------------------------------------------------------------------------

library(DropletUtils)
library(Seurat)
source(file = snakemake@params[["sce_functions"]])

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_10"]])
separate_fractions_clustering <- snakemake@params[["separate_fractions_clustering"]]

#-------------------------------------------------------------------------------
# Clustering

if(separate_fractions_clustering){
  print("separated fractions")
  
  # separate fractions for separate clustering 
  sce_hsc <- sce[,sce$Fraction_ID == "hsc"]
  sce_str <- sce[,sce$Fraction_ID == "str"]
  
  # convert so seurat objects
  if(colData(sce)$Correction_method[1] == "Seurat"){ # if seurat was used for BC
    seurat_hsc <- as.Seurat(sce_hsc, counts = "counts", data = "corrected")
    seurat_str <- as.Seurat(sce_str, counts = "counts", data = "corrected")
  }else if(colData(sce)$Correction_method[1] == "FastMNN"){ # if mnncorrect was used
    seurat_hsc <- as.Seurat(sce_hsc, counts = "counts", data = "reconstructed")
    seurat_str <- as.Seurat(sce_str, counts = "counts", data = "reconstructed")
  }
  
  # clustering based on PCA on batch corrected values
  seurat_hsc <- FindNeighbors(seurat_hsc, dims = 1:10, reduction = "PCA")
  seurat_hsc <- FindClusters(seurat_hsc, resolution = 0.5)

  seurat_str <- FindNeighbors(seurat_str, dims = 1:10, reduction = "PCA")
  seurat_str <- FindClusters(seurat_str, resolution = 0.5)

  # transfer labels back and add together
  sce_hsc$cluster_seurat <- Idents(seurat_hsc)
  sce_hsc$cluster_seurat <- as.numeric(unfactor(sce_hsc$cluster_seurat))
  sce_str$cluster_seurat <- Idents(seurat_str)
  sce_str$cluster_seurat <- as.numeric(sce_str$cluster_seurat) + 
    length(unique(sce_hsc$cluster_seurat))
  sce <- cbind(sce_hsc, sce_str)
  sce$cluster_seurat <- factor(sce$cluster_seurat, 
                               levels = sort(unique(sce$cluster_seurat)))
  sce$seurat_mode <- rep("fractions_separated", ncol(sce))
  
}else{
  
  if(sce$Correction_method[1] == "Seurat"){ # if seurat was used for BC
    seurat <- as.Seurat(sce, counts = "counts", data = "corrected")
  }else if(sce$Correction_method[1] == "FastMNN"){ # if mnncorrect was used
    seurat <- as.Seurat(sce, counts = "counts", data = "reconstructed")
  }

  seurat <- FindNeighbors(seurat, dims = 1:10, reduction = "PCA")
  seurat <- FindClusters(seurat, resolution = 0.5)
  sce$cluster_seurat <- Idents(seurat)
  sce$seurat_mode <- rep("fractions_together", ncol(sce))
  
}
 
print(sce)
saveRDS(sce, file = snakemake@output[["sce_11"]])