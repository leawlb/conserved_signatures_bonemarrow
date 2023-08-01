#-------------------------------------------------------------------------------

library(SingleCellExperiment)
library(mclust, quietly = TRUE) 

set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_input"]])

cluster_curr <- snakemake@wildcards[["cluster"]]
fraction_curr <- snakemake@wildcards[["fraction"]]

if(cluster_curr == 2){
  sce_clust <- sce[,sce$cluster_louvain %in% c(2, 4)]
  cluster_curr <- "2_4"
}else{
  sce_clust <- sce[,sce$cluster_louvain == cluster_curr]
}

print(cluster_curr)
sce_clust$subcluster <- vector(length = ncol(sce_clust))
sce_clust$subcluster_annotation <- vector(length = ncol(sce_clust))

#-------------------------------------------------------------------------------

genes_list_shared <- readRDS(file = snakemake@input[["genes_list_shared"]])

celltypes <- unique(sce_clust$annotation_cluster)
print(celltypes)

if(length(celltypes) == 1){
  genes_shared <- genes_list_shared[[celltypes]]$at_least_three
}else if(length(celltypes) == 2){
  genes_shared <- intersect(genes_list_shared[[celltypes[1]]]$at_least_three, 
                            genes_list_shared[[celltypes[2]]]$at_least_three)
}

#-------------------------------------------------------------------------------

genes_subcl <- read.csv(file = snakemake@input[["gene_list_subcl"]], 
                        header = TRUE, 
                        sep = ";", 
                        check.names=FALSE, 
                        stringsAsFactors=FALSE, 
                        as.is=TRUE, 
                        colClasses = "character")

genes_subcl <- genes_subcl[genes_subcl$fraction == fraction_curr,]
genes_subcl <- genes_subcl[genes_subcl$cluster == cluster_curr,]

subcl <- unique(genes_subcl$gene[genes_subcl$purpose == "subclustering"])

stopifnot(subcl %in% genes_shared)

print(paste("subcl genes used:", subcl[which(subcl %in% genes_shared)]))
print(paste("subcl genes not shared, not used:", 
            subcl[which(!subcl %in% genes_shared)]))

genes_all <- unique(genes_subcl$gene)
poss_subcl <- genes_all[which(genes_all %in% genes_shared)]

print(paste("subcl genes possible:", poss_subcl))

#-------------------------------------------------------------------------------

anno_subcl <- read.csv(file = snakemake@input[["anno_subcl"]], 
                       header = TRUE, 
                       sep = ";", 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       as.is=TRUE, 
                       colClasses = "character")

anno_subcl <- anno_subcl[anno_subcl$fraction == fraction_curr,]
anno_subcl <- anno_subcl[anno_subcl$cluster == cluster_curr,]

nr_subcl <- max(anno_subcl$subcluster)
print(cluster_curr)
print(nr_subcl)

#-------------------------------------------------------------------------------

data <- logcounts(sce_clust)
print(paste("final print", subcl))
data <- t(data[rownames(data) %in% subcl,])
head(data)
clust <- Mclust(data, G = nr_subcl,  modelNames = c("EII",
                                                    "VII",
                                                    "EEI",
                                                    "VEI",
                                                    "EVI",
                                                    "VVI",
                                                    "VEE",
                                                    "EVE",
                                                    "VVE",
                                                    "EEV",
                                                    "VEV",
                                                    "EVV"))
# all multivariate models for data with observations > variables, except VVV
# VVV doesn't work on the cluster, but locally run models all contained one E 

summary(clust)

sce_clust$subcluster <- as.character(clust$classification)
print(unique(sce_clust$subcluster))
print(anno_subcl)
 
for(i in 1:nr_subcl){
  sce_clust$annotation_subcluster[
    sce_clust$subcluster == i] <- anno_subcl$annotation_subcluster[
      anno_subcl$subcluster == i]
}
print(unique(sce_clust$annotation_subcluster))
#-------------------------------------------------------------------------------

saveRDS(sce_clust, snakemake@output[["sce_output"]])