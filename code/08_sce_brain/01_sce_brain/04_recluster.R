library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
set.seed(37)

#base_path_temp <- snakemake@params[["base_path"]] # try and see what happens if base_path is used --- I don't have access to OE0433
base_path <- snakemake@params[["base_path"]]
brain_path <- snakemake@params[["brain_path"]]

data <- readRDS(snakemake@input[["data_input"]])

data.updated <- UpdateSeuratObject(object = data)  # available data is v3 Seurat

nDEGs <- readRDS(paste0(base_path, brain_path, "01_list_nDEGs_all.rds"))
cons_markers <- readRDS(paste0(base_path, brain_path, "02_marker_conserved_primates.rds"))

# clustering resolutions are based on values optimized using conserved markers
resolution <- data.frame(species = c("human", "marmoset", "macaque", "mouse"),
                         resolution = c(0.5, 0.6, 0.7, 0.7)) # chosen in report

# table(data.updated$subclass_label, data.updated$orig.ident)
#            human macaque marmoset mouse
# L2/3 IT     1583       0     1600  1599
# L5 ET        505     483      499   800
# L5 IT       2347    1984     2400  2282
# L5/6 NP      825     435      850   983
# L6 CT       1110     893     1200   912
# L6 IT        600     600      600   600
# L6 IT Car3   200     200      200    13
# L6b          635     242      817  1053


conserved_signature <- lapply(names(cons_markers), function(x){
  markers <- cons_markers[[x]]
  ndeg <- nDEGs[[x]]
  return(markers[markers %in% ndeg])
})
names(conserved_signature) <- names(cons_markers)


for(sp in c("human", "marmoset", "macaque", "mouse")){
  species <- subset(data.updated, subset = orig.ident == sp)
  
  DefaultAssay(species) <- "RNA"
  species <- ScaleData(species)
  species <- RunPCA(species,
                  features = unique(unlist(conserved_signature)))
  species <- RunUMAP(species,
                   reduction = "pca",
                   dims = 1:50,
                   n.components = 2)
  species <- FindNeighbors(species,
                         reduction = "pca",
                         dims = 1:50)
  species <- FindClusters(species, 
                          resolution = resolution[which(resolution$species == sp), "resolution"])
  
  species@commands <- list()
  
  saveRDS(species,
          file = paste0(base_path, brain_path, "04_recluster/",
                        sp, "_reclust_sig.rds"))
  
  species <- RunPCA(species,
                    features = unique(unlist(cons_markers)))
  species <- RunUMAP(species,
                     reduction = "pca",
                     dims = 1:50,
                     n.components = 2)
  species <- FindNeighbors(species,
                           reduction = "pca",
                           dims = 1:50)
  species <- FindClusters(species, 
                          resolution = resolution[which(resolution$species == sp), "resolution"])
  
  species@commands <- list()
  
  saveRDS(species,
          file = paste0(base_path, brain_path, "04_recluster/",
                        sp, "_reclust_core.rds"))
}





# identify proper order of patterns for final visualization
for(sp in c("human", "marmoset", "macaque", "mouse")){
  for(genes in c("core", "sig")){
    data <- readRDS(paste0(base_path, brain_path, "04_recluster/",
                           sp, "_reclust_", genes, ".rds"))
    
    plot <- data@meta.data[,c("subclass_label", "seurat_clusters")] %>%
      table() %>% as.matrix()
    
    plot_gg <- plot %>%
      as.data.frame() %>%
      group_by(seurat_clusters) %>%
      mutate(clust_prop = Freq/sum(Freq)) %>%
      ungroup %>%
      group_by(subclass_label) %>%
      mutate(subclass_prop = Freq/sum(Freq))
    
    print(
      ggplot(plot_gg, 
             aes(x = factor(subclass_label, 
                            levels = c("L2/3 IT", "L5 IT", "L6 IT", "L6 IT Car3", 
                                       "L5 ET", "L5/6 NP", "L6 CT", "L6b")), 
                 y = seurat_clusters, 
                 fill= clust_prop*100)) + 
        geom_tile() +
        scale_fill_continuous("% cells/cluster", 
                              limits=c(0, 100), 
                              breaks=seq(0,100,by=20),
                              low = "white", high = "blue") +
        theme_classic()+
        theme(axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
        labs(y  = "new clusters", x = "original cell types", title = paste(sp, genes))
    )
  }
}





############################################
## cluster mouse using only human markers ##
############################################
markers_conservation <- readRDS(paste0(base_path, brain_path, "02_marker_expression_primates.rds"))
human_markers <- do.call(rbind, markers_conservation) %>%
  as.data.frame() %>%
  filter(!is.na(human)) %>% 
  rownames() %>%
  sub("\\..*", "", .) %>%
  unique()

species <- subset(data.updated, subset = orig.ident == "mouse")

DefaultAssay(species) <- "RNA"
species <- ScaleData(species)
species <- RunPCA(species,
                  features = human_markers)
species <- RunUMAP(species,
                   reduction = "pca",
                   dims = 1:50,
                   n.components = 2)
species <- FindNeighbors(species,
                         reduction = "pca",
                         dims = 1:50)
species <- FindClusters(species, 
                        resolution = resolution[which(resolution$species == "mouse"), "resolution"])

species@commands <- list()

saveRDS(species, snakemake@output[["species"]])
