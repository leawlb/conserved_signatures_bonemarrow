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

# assign cluster annotations from cluster_annotations.txt
sce$annotation_cluster <- vector(length = ncol(sce))

for(i in unique(clust_anno$Cluster)){
  sce$annotation_cluster[sce$cluster_louvain == i] <- clust_anno$Annotation[clust_anno$Cluster == i]
}

sce <- sce[,-which(sce$annotation_cluster == "remove")]

sce$annotation_cluster <- factor(sce$annotation_cluster,
                                 levels = clust_anno$Annotation[clust_anno$Annotation != "remove"])

saveRDS(sce, file = snakemake@output[["sce_output"]])
