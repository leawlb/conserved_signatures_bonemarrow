# visualize expression of conserved signature genes along pseudotime
# the content is now in cellrank_report_sep for automatic visualisation with snakemake

library(Seurat)
library(tidyverse)

markers <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own/01_gens/geneset_list_hsc")

# load all three trajectories
ery_pseudo <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/05_cellrank/05_bsce/sce_ery")
ery_pseudo <- CreateSeuratObject(counts = ery_pseudo@assays@data@listData[["logcounts"]],
                                 meta.data = ery_pseudo@colData@listData)

lym_pseudo <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/05_cellrank/05_bsce/sce_lym")
lym_pseudo <- CreateSeuratObject(counts = lym_pseudo@assays@data@listData[["logcounts"]],
                                 meta.data = lym_pseudo@colData@listData)

neu_pseudo <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/05_cellrank/05_bsce/sce_neu")
neu_pseudo <- CreateSeuratObject(counts = neu_pseudo@assays@data@listData[["logcounts"]],
                                 meta.data = neu_pseudo@colData@listData)

# demonstrates that each species has approximate same distribution of
# cell types across pseudotime (at least within largest, ery, trajectory)
ggplot(ery_pseudo@meta.data,
       aes( #x = pseudotime,
         x = rank(pseudotime),
         fill = celltypes,
         color = celltypes)) +
  facet_wrap(~Species_ID,
             ncol = 1,
             scales = "free_y") +
  geom_density(alpha = 0.4, linewidth = 0.2) +
  theme_classic() +
  labs(y = NULL)


plot_data <- ery_pseudo@meta.data %>%
  select("pseudotime", "celltypes", "Mouse_ID", 
         "Species_ID", "Age_weeks", "Age_ID")
for(n in 1:length(markers)){
  cell_sp_markers <- markers[[n]][["conserved_signature"]]
  for(g in cell_sp_markers){
    plot_data[,g] <-unlist(FetchData(ery_pseudo, g))
  }
}
plot_data[,7:ncol(plot_data)] <- scale(plot_data[,7:ncol(plot_data)])

pdf("/omics/odcf/analysis/OE0538_projects/DO-0008/data/test_reproducibility4/main_analysis/sce_objects/03_sce_analysis/05_cellrank/06_pseudotime_figs/marker_genes_pseudo_ery.pdf",
    width = 5, height = 4)
for(n in 1:length(markers)){
  cell_sp <- markers[[n]][["conserved_signature"]]
  plot_genes <- plot_data[,c(colnames(plot_data)[1:6],cell_sp)]
  plot_genes <- gather(plot_genes, gene, expression, 
                       -pseudotime, -celltypes, -Mouse_ID, 
                       -Species_ID, -Age_weeks, -Age_ID)
  print(
    ggplot(plot_genes,
    aes(x = rank(pseudotime),
        y = expression,
        color = gene)) +
    geom_smooth(se = F) +
    guides(color = "none") +
    labs(x = "Pseudotime Ery", y = "Scaled Expression",
         title = names(markers)[n]) +
    theme_classic() +
    theme(axis.text = element_blank())
  )
}
dev.off()



plot_data <- lym_pseudo@meta.data %>%
  select("pseudotime", "celltypes", "Mouse_ID", 
         "Species_ID", "Age_weeks", "Age_ID")
for(n in 1:length(markers)){
  cell_sp_markers <- markers[[n]][["conserved_signature"]]
  for(g in cell_sp_markers){
    plot_data[,g] <-unlist(FetchData(lym_pseudo, g))
  }
}
plot_data[,7:ncol(plot_data)] <- scale(plot_data[,7:ncol(plot_data)])

pdf("/omics/odcf/analysis/OE0538_projects/DO-0008/data/test_reproducibility4/main_analysis/sce_objects/03_sce_analysis/05_cellrank/06_pseudotime_figs/marker_genes_pseudo_lym.pdf",
    width = 5, height = 4)
for(n in 1:length(markers)){
  cell_sp <- markers[[n]][["conserved_signature"]]
  plot_genes <- plot_data[,c(colnames(plot_data)[1:6],cell_sp)]
  plot_genes <- gather(plot_genes, gene, expression, 
                       -pseudotime, -celltypes, -Mouse_ID, 
                       -Species_ID, -Age_weeks, -Age_ID)
  print(
    ggplot(plot_genes,
           aes(x = rank(pseudotime),
               y = expression,
               color = gene)) +
      geom_smooth(se = F) +
      guides(color = "none") +
      labs(x = "Pseudotime Lym", y = "Scaled Expression",
           title = names(markers)[n]) +
      theme_classic() +
      theme(axis.text = element_blank())
  )
}
dev.off()



plot_data <- neu_pseudo@meta.data %>%
  select("pseudotime", "celltypes", "Mouse_ID", 
         "Species_ID", "Age_weeks", "Age_ID")
for(n in 1:length(markers)){
  cell_sp_markers <- markers[[n]][["conserved_signature"]]
  for(g in cell_sp_markers){
    plot_data[,g] <-unlist(FetchData(neu_pseudo, g))
  }
}
plot_data[,7:ncol(plot_data)] <- scale(plot_data[,7:ncol(plot_data)])

pdf("/omics/odcf/analysis/OE0538_projects/DO-0008/data/test_reproducibility4/main_analysis/sce_objects/03_sce_analysis/05_cellrank/06_pseudotime_figs/marker_genes_pseudo_neu.pdf",
    width = 5, height = 4)
for(n in 1:length(markers)){
  cell_sp <- markers[[n]][["conserved_signature"]]
  plot_genes <- plot_data[,c(colnames(plot_data)[1:6],cell_sp)]
  plot_genes <- gather(plot_genes, gene, expression, 
                       -pseudotime, -celltypes, -Mouse_ID, 
                       -Species_ID, -Age_weeks, -Age_ID)
  print(
    ggplot(plot_genes,
           aes(x = rank(pseudotime),
               y = expression,
               color = gene)) +
      geom_smooth(se = F) +
      guides(color = "none") +
      labs(x = "Pseudotime Neu", y = "Scaled Expression",
           title = names(markers)[n]) +
      theme_classic() +
      theme(axis.text = element_blank())
  )
}
dev.off()


pdf("/omics/odcf/analysis/OE0538_projects/DO-0008/data/test_reproducibility4/main_analysis/sce_objects/03_sce_analysis/05_cellrank/06_pseudotime_figs/marker_genes_ery_species_pseudotime.pdf",
    width = 7, height = 4)
for(g in sort(unique(unlist(lapply(markers, 
                                    function(x){x["conserved_signature"]}))))) {
  plot_data <- ery_pseudo@meta.data %>%
    select("pseudotime", "celltypes", "Mouse_ID", 
           "Species_ID", "Age_weeks", "Age_ID")
  plot_data$gene <-unlist(FetchData(ery_pseudo, g))
  
  print(
    ggplot(plot_data,
           aes(x = rank(pseudotime),
               y = gene,
               color = Species_ID)) +
      geom_smooth() +
      theme_classic() +
      labs(x = "Pseudotime Ery", y = g)  +
      theme(axis.text = element_blank())
  )
}
dev.off()


pdf("/omics/odcf/analysis/OE0538_projects/DO-0008/data/test_reproducibility4/main_analysis/sce_objects/03_sce_analysis/05_cellrank/06_pseudotime_figs/marker_genes_lym_species_pseudotime.pdf",
    width = 7, height = 4)
for(g in sort(unique(unlist(lapply(markers, 
                                   function(x){x["conserved_signature"]}))))) {
  plot_data <- lym_pseudo@meta.data %>%
    select("pseudotime", "celltypes", "Mouse_ID", 
           "Species_ID", "Age_weeks", "Age_ID")
  plot_data$gene <-unlist(FetchData(lym_pseudo, g))
  
  print(
    ggplot(plot_data,
           aes(x = rank(pseudotime),
               y = gene,
               color = Species_ID)) +
      geom_smooth() +
      theme_classic() +
      labs(x = "Pseudotime Lym", y = g)  +
      theme(axis.text = element_blank())
  )
}
dev.off()


pdf("/omics/odcf/analysis/OE0538_projects/DO-0008/data/test_reproducibility4/main_analysis/sce_objects/03_sce_analysis/05_cellrank/06_pseudotime_figs/marker_genes_neu_species_pseudotime.pdf",
    width = 7, height = 4)
for(g in sort(unique(unlist(lapply(markers, 
                                   function(x){x["conserved_signature"]}))))) {
  plot_data <- neu_pseudo@meta.data %>%
    select("pseudotime", "celltypes", "Mouse_ID", 
           "Species_ID", "Age_weeks", "Age_ID")
  plot_data$gene <-unlist(FetchData(neu_pseudo, g))
  
  print(
    ggplot(plot_data,
           aes(x = rank(pseudotime),
               y = gene,
               color = Species_ID)) +
      geom_smooth() +
      theme_classic() +
      labs(x = "Pseudotime Neu", y = g)  +
      theme(axis.text = element_blank())
  )
}
dev.off()
