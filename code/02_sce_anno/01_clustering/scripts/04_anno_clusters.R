#-------------------------------------------------------------------------------
library(DropletUtils, quietly = TRUE) 

#TODO: either completely delete when not needed anymore or take from CSV

sce <- readRDS(file = snakemake@input[["sce_input"]])

anno_clusters <- read.csv(file = snakemake@input[["anno_clusters"]],
                       header = TRUE, 
                       sep = ";", 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       as.is=TRUE, 
                       colClasses = "character")

fraction_curr <- sce$Fraction_ID[1]
anno_clusters <- anno_clusters[anno_clusters$Fraction == fraction_curr,]

# assign cluster annotations from cluster_annotations.txt
sce$annotation_cluster <- vector(length = ncol(sce))

for(i in unique(anno_clusters$Cluster)){
  sce$annotation_cluster[sce$cluster_louvain == i] <- anno_clusters$Annotation[anno_clusters$Cluster == i]
}

sce <- sce[,-which(sce$annotation_cluster == "remove")]

sce$annotation_cluster <- factor(sce$annotation_cluster,
                                 levels = anno_clusters$Annotation[anno_clusters$Annotation != "remove"])

table(is.na(sce$annotation_cluster))

stopifnot(!is.na(sce$annotation_cluster))
stopifnot(!is.na(sce$cluster_louvain))

saveRDS(sce, file = snakemake@output[["sce_output"]])
