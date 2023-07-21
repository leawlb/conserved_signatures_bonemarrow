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

genes_shared_list <- readRDS(file = snakemake@input[["genes_shared_list"]])

celltypes <- unique(sce_clust$annotation_cluster)
print(celltypes)

if(length(celltypes) == 1){
  genes_shared <- genes_shared_list[[celltypes]]$at_least_three
}else if(length(celltypes) == 2){
  genes_shared <- intersect(genes_shared_list[[celltypes[1]]]$at_least_three, 
                            genes_shared_list[[celltypes[2]]]$at_least_three)
}

#-------------------------------------------------------------------------------

subcl_genes <- read.csv(file = snakemake@input[["subclustering_genes"]], 
                        header = TRUE, 
                        sep = ";", 
                        check.names=FALSE, 
                        stringsAsFactors=FALSE, 
                        as.is=TRUE, 
                        colClasses = "character")

subcl_genes <- subcl_genes[subcl_genes$fraction == fraction_curr,]
subcl_genes <- subcl_genes[subcl_genes$cluster == cluster_curr,]

subcl <- unique(subcl_genes$gene[subcl_genes$purpose == "subclustering"])

stopifnot(subcl %in% genes_shared)
print(subcl)

#-------------------------------------------------------------------------------

subcl_anno <- read.csv(file = snakemake@input[["subclustering_annotation"]], 
                       header = TRUE, 
                       sep = ";", 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       as.is=TRUE, 
                       colClasses = "character")

subcl_anno <- subcl_anno[subcl_anno$fraction == fraction_curr,]
subcl_anno <- subcl_anno[subcl_anno$cluster == cluster_curr,]

nr_subcl <- max(subcl_anno$subcluster)
print(cluster_curr)
print(nr_subcl)

#-------------------------------------------------------------------------------

data <- logcounts(sce_clust)
data <- t(data[rownames(data) %in% subcl,])
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
print(subcl_anno)
 
for(i in 1:nr_subcl){
  sce_clust$annotation_subcluster[
    sce_clust$subcluster == i] <- subcl_anno$annotation_subcluster[
      subcl_anno$subcluster == i]
}
print(unique(sce_clust$annotation_subcluster))
#-------------------------------------------------------------------------------

saveRDS(sce_clust, snakemake@output[["sce_output"]])