#-------------------------------------------------------------------------------

library(monocle3)
library(scran)

set.seed(37)

#-------------------------------------------------------------------------------
sce <- readRDS(snakemake@input[["sce_input"]])

# cell cycle genes taken from go term cell cycle process, ID = GO:0022402
# downloaded from https://www.informatics.jax.org/go/term/GO:0022402 on 2024-01-12
cell_cycle_df <- read.csv(file = snakemake@input[["cell_cycle_genes"]], sep = ",")

#sce <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10")

#-------------------------------------------------------------------------------
# prepare SCE object for conversion

cell_cycle_genes <- unique(cell_cycle_df$Symbol)
print(length(cell_cycle_genes))

sce <- sce[which(!rownames(sce) %in% cell_cycle_genes),]

colData(sce) <- colData(sce)[,colnames(colData(sce)) %in% c("individual", 
                                                            "Species_ID",
                                                            "Age_ID",
                                                            "Fraction_ID",
                                                            "celltypes", 
                                                            "category")]
rowData(sce) <- rowData(sce)[,colnames(rowData(sce)) %in% c("ID", 
                                                            "Symbol")]
rowData(sce)$gene_short_name <- rowData(sce)$Symbol


#-------------------------------------------------------------------------------
# conversion

expression_matrix <- counts(sce)
cell_metadata <- colData(sce)
gene_annotation <- rowData(sce)
  
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

print("conversion done")

#-------------------------------------------------------------------------------
# main steps for monocle3 prep, using default options

# pre-processing
cds <- preprocess_cds(cds, num_dim = 50)
  
# batch correction (on individual, using MNN)
cds <- align_cds(cds, alignment_group = "individual")

cds_list <- list()
cds_list[["Mk/Ery"]] <- cds[,cds$category %in% c("Stem", "MPPs", "Cycling", "Mk/Ery")]
cds_list[["Myeloid"]] <- cds[,cds$category %in% c("Stem", "MPPs", "Cycling", "Myeloid")]
cds_list[["Lymphoid"]] <- cds[,cds$category %in% c("Stem", "MPPs", "Cycling", "Lymphoid")]
  
cds_list <- lapply(cds_list, function(cds){
  
# dimensionality reduction (3D)
  cds <- reduce_dimension(cds, max_components = 2)

#reducedDims(cds)$UMAP <- reducedDims(sce)$UMAP
#reducedDims(cds)$PCA <- reducedDims(sce)$PCA

# cluster cells
# 0.0003
  cds <- cluster_cells(cds, cluster_method = "leiden", resolution = 0.004,
                       random_seed = 373, reduction_method = "UMAP")
  print(cds)
  return(cds)
})
# save image
#pdf(snakemake@output[["pdf_output1"]])
#plot_cells_3d(cds, label_principal_points = TRUE, color_cells_by = "partition")
#dev.off() 

# save image 
#pdf(snakemake@output[["pdf_output2"]])
#plot_cells_3d(cds, label_principal_points = TRUE, color_cells_by = "celltypes")
#dev.off() 

saveRDS(cds_list, snakemake@output[["cds_output"]])
