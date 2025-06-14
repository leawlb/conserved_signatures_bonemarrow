library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
set.seed(37)

data <- readRDS(snakemake@input[["data_input"]])
data.updated <- UpdateSeuratObject(object = data)  # available data is v3 Seurat

nDEGs <- readRDS(snakemake@input[["nDEGS"]])
cons_markers <- readRDS(snakemake@input[["cons_markers"]])

base_path <- snakemake@params[["base_path"]]
brain_path <- snakemake@params[["brain_path"]]

# clustering resolutions are based on values optimized using conserved markers
# reclustering with conserved markers not currently required
# resolution_marker <- data.frame(
#   species = c("human", "marmoset", "macaque", "mouse"),
#   resolution = c(0.5, 0.6, 0.7, 0.7)) # chosen in report

# and for signature
resolution_sign <- data.frame(
  species = c("human", "marmoset", "macaque", "mouse"),
  resolution = c(1, 1.1, 0.8, 1)) # chosen in report


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
  
  print(sp)
  
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
                          resolution = resolution_sign[
                            which(resolution_sign$species == sp), "resolution"])
  
  species@commands <- list()
  
  saveRDS(species,
          file = paste0(base_path, brain_path, "04_recluster/",
                        sp, "_reclust_sig.rds"))
  
  # reclustering with conserved markers not currently required
  # species <- RunPCA(species,
  #                   features = unique(unlist(cons_markers)))
  # species <- RunUMAP(species,
  #                    reduction = "pca",
  #                    dims = 1:50,
  #                    n.components = 2)
  # species <- FindNeighbors(species,
  #                          reduction = "pca",
  #                          dims = 1:50)
  # species <- FindClusters(species, 
  #                         resolution = resolution_marker[
  #                           which(resolution_marker$species == sp), "resolution"])
  # 
  # species@commands <- list()
  # 
  # saveRDS(species,
  #         file = paste0(base_path, brain_path, "04_recluster/",
  #                       sp, "_reclust_core.rds"))
}

print("signature reclustering done")

############################################
## cluster mouse using only human markers ##
############################################

markers_conservation <- readRDS(snakemake@input[["markers_conservation"]])
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
                        resolution = resolution_sign[
                          which(resolution_sign$species == "mouse"),
                          "resolution"])

species@commands <- list()

saveRDS(species, snakemake@output[["species_human_only"]])
