#-------------------------------------------------------------------------------
library(DropletUtils, quietly = TRUE) 

#TODO: either completely delete when not needed anymore or take from CSV

sce <- readRDS(file = snakemake@input[["sce_input"]])

clust_anno <- read.csv(file = snakemake@input[["cluster_annotations"]],
                       header = TRUE, 
                       sep = ";", 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       as.is=TRUE, 
                       colClasses = "character")

fraction_curr <- sce$Fraction_ID[1]
clust_anno <- clust_anno[clust_anno$Fractions == "fraction_curr",]

# assign cluster annotations
sce$annotation_cluster <- vector(length = ncol(sce))

for(i in unique(clust_anno$Cluster)){
  sce$annotation_cluster[sce$cluster_louvain == i] <- clust_anno$Annotation[clust_anno$Cluster == i]
}

if(fraction_curr == "hsc"){
  
  sce$annotation_cluster <- factor(sce$annotation_cluster,
                                   levels = c(
                                     "HSCs",
                                     "MPPs early",
                                     "MPPs differentiating",
                                     "Myeloid/Granu early",
                                     "Myeloid/Granu late",
                                     "Baso/Eo/Mast progs.",
                                     "Mk",
                                     "Mk/Erythroid progs.",
                                     "Erythroid",
                                     "Cycling1",
                                     "Cycling2"
                                     ))

}else if(fraction_curr == "str"){

  
  sce$annotation_cluster <- factor(sce$annotation_cluster,
                                 levels = c(
                                   "Mesenchymal",
                                   "Fibroblasts",
                                   "CAR cells 1",
                                   "CAR cells 2",
                                   "Arteriolar/Capillary ECs",
                                   "Arteriolar/Mixed ECs",
                                   "Sinusoidal ECs",
                                   "Pericytes/Smooth muscle cells",
                                   "Skeletal muscle cells",
                                   "mono lineage 1",
                                   "mono lineage 2",
                                   "neutro lineage",
                                   "lymphoid 1",
                                   "lymphoid 2",
                                   "lymphoid 3",
                                   "ery lineage",
                                   "mast/eo or granu"
                                 ))
}

saveRDS(sce, file = snakemake@output[["sce_output"]])
