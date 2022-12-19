#-------------------------------------------------------------------------------

library(DropletUtils)
library(Seurat)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_10"]])
remove_mature_cells_clustering <- snakemake@params[["remove_mature_cells_clustering"]]
separate_fractions_clustering <- snakemake@params[["separate_fractions_clustering"]]

# prepare SCE
sce <- get_mature_cts_all(sce)
if(remove_mature_cells_clustering){
  sce <- sce[,sce$mature_cells == "FALSE"]
}

# separate fractions for separate clustering 
sce_hsc <- sce[,sce$Fraction_ID == "hsc"]
sce_str <- sce[,sce$Fraction_ID == "str"]

# find contamination but do not remove it (save it for later)
sce_hsc <- find_contamination_ref_all(sce_hsc, input = "hsc")
sce_str <- find_contamination_ref_all(sce_str, input = "stromal")

#-------------------------------------------------------------------------------
# Clustering

if(separate_fractions_clustering){
  print("separated fractions")
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
  sce <- cbind(sce_hsc, sce_str)
  if(bc_use == "seurat3"){ # if seurat was used for BC
    seurat <- as.Seurat(sce, counts = "counts", data = "corrected")
  }else if(bc_use == "mnncorrect"){ # if mnncorrect was used
    seurat <- as.Seurat(sce, counts = "counts", data = "reconstructed")
  }

  seurat <- FindNeighbors(seurat, dims = 1:10, reduction = "PCA")
  seurat <- FindClusters(seurat, resolution = 0.5)
  sce$cluster_seurat <- Idents(seurat)
  sce$seurat_mode <- rep("fractions_together", ncol(sce))
  
}
 
print(sce)
saveRDS(sce, file = snakemake@output[["sce_11"]])