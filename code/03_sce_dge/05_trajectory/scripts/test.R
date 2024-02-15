

library(tidyverse)
library(monocle3)
library(scran)



expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))


### basic workflow with example DS

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "batch")

cds <- reduce_dimension(cds)
plot_cells(cds)



cds <- cluster_cells(cds)

cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
cds
reducedDims(cds)$UMAP
reducedDims(cds)$PCA

#-------------------------------------------------------------------------------

sce <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/02_sce_anno/10_anns/sce_str-10")
sce <- sce[,which(sce$Species_ID == "mmus" & sce$Age_ID == "yng")]
colData(sce) <- colData(sce)[,colnames(colData(sce)) %in% c("individual", 
                                                            "Species_ID",
                                                            "Age_ID",
                                                            "Fraction_ID",
                                                            "celltypes")]
rowData(sce) <- rowData(sce)[,colnames(rowData(sce)) %in% c("ID", 
                                                            "Symbol")]
rowData(sce)$gene_short_name <- rowData(sce)$Symbol

expression_matrix <- counts(sce)
cell_metadata <- colData(sce)
gene_annotation <- rowData(sce)

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "individual")

cds <- reduce_dimension(cds)

plot_cells(cds)


cds <- cluster_cells(cds)

cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds,  label_principal_points = TRUE)

pt <- pseudotime(cds, reduction_method = "UMAP")


reducedDims(cds)$Aligned
sce_hsc$pseudotime <- pt

source("/home/l012t/repositories/Interspecies_BM_phd/code/source/plotting.R")

library(viridis)

umap_base(sce_hsc, "pseudotime")+
  scale_color_viridis("")

median(sce_hsc[,sce_hsc$celltypes == "Monocyte progs."]$pseudotime)
mean(sce_hsc[,sce_hsc$celltypes == "Mk progs."]$pseudotime)


sce_hsc$category

sce_hsc <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10")
sce_str <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/02_sce_anno/10_anns/sce_str-10")


sort(unique(sce_str$Object_ID))


cds <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/03_sce_analysis/11_mprp/cds_hsc")
cds <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/03_sce_analysis/12_mgrp/cds_hsc")
cds <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/03_sce_analysis/13_pstm/cds_hsc")

colData(cds_hsc)

core_cons_list_hsc <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/03_sce_analysis/07_core/core_cons_list_hsc")
core_cons_list_hsc$HSCs

core_cons_list_hsc_new <- list()
for(i in 1:length(core_cons_list_hsc)){
  core_cons_list_hsc_new[[i]] <- core_cons_list_hsc[[i]]$marker_cons_ct
}

saveRDS(core_cons_list_hsc_new, "/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/03_sce_analysis/07_core/core_cons_list_hsc_manual")


library("xlsx")


cell_cycle_df <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/03_sce_analysis/cell_cycle_genes/go_term_cellcycleprocess_mouse.csv", sep = ",")

cell_cycle_genes <- unique(cell_cycle_df$Symbol)




plot_cells_3d(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE)

plot_cells(cds, color_cells_by = "celltypes")




print(length(cell_cycle_genes))





colors_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/colors/colors.txt"
source("/home/l012t/repositories/Interspecies_BM_phd/code/source/colors.R")

umap_base(sce_hsc[,sample(c(1:ncol(sce_hsc)), ncol(sce_hsc))], "Species_ID")+
  scale_color_manual(values = col_spc)

