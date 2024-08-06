#-------------------------------------------------------------------------------
# sub-cluster selected clusters using mclust and a defined set of genes

library(mclust, quietly = TRUE) 
set.seed(37)

#-------------------------------------------------------------------------------
# load and subset sce object

sce <- base::readRDS(file = snakemake@input[["sce_input"]])
sep <- base::readRDS(file = snakemake@input[["sep"]]) # only the cluster nr

anno_subcl_path <- snakemake@input[["anno_subcl"]]
gene_list_shared_path <- snakemake@input[["genes_list_shared"]]
gene_list_subcl_path <- snakemake@input[["gene_list_subcl"]]

cluster_curr <- snakemake@wildcards[["cluster"]]
fraction_curr <- snakemake@wildcards[["fraction"]]

print(sce)
print(cluster_curr)

#-------------------------------------------------------------------------------

# hsc clusters 2 and are subclustered together
# this script only runs for cluster "2", but cluster 4 is included below (HSPCs)
if(cluster_curr == 2){
  sce_clust <- sce[,which(sce$cluster_louvain %in% c(2, 4))]
  cluster_curr <- "2_4"
}else{
  sce_clust <- sce[,which(sce$cluster_louvain == cluster_curr)]
}

sce_clust$subcluster <- vector(length = ncol(sce_clust))

#-------------------------------------------------------------------------------
# load and subset shared nDGEs (shared between three of four species)

genelist_shared <- base::readRDS(file = gene_list_shared_path)

celltypes <- base::unique(sce_clust$annotation_cluster)
print(celltypes)

if(length(celltypes) == 1){
  genes_shared <- genelist_shared[[celltypes]]$three
  
  # for subclustering clusters 2 and 4, use intersection of both clusters
}else if(length(celltypes) == 2){ 
  genes_shared <- base::intersect(genelist_shared[[celltypes[1]]]$three, 
                                  genelist_shared[[celltypes[2]]]$three)
}

#-------------------------------------------------------------------------------
# load and subset genes for sub-clustering (picked out manually)
# cluster is currently 2_4

genes_subcl <- utils::read.csv(file = gene_list_subcl_path, 
                               header = TRUE, 
                               sep = ";", 
                               check.names=FALSE, 
                               stringsAsFactors=FALSE, 
                               as.is=TRUE, 
                               colClasses = "character")
print(head(genes_subcl))

genes_subcl <- genes_subcl[genes_subcl$fraction == fraction_curr,]
genes_subcl <- genes_subcl[genes_subcl$cluster == cluster_curr,]

subcl <- base::unique(genes_subcl$gene[genes_subcl$purpose == "subclustering"])

#-------------------------------------------------------------------------------
# test that all subclustering genes are also shared nDGEs (failsave)
# shared between three of four species as determined in 02_nDGE

print(base::paste("subcl genes not shared, not used:", 
                  subcl[which(!subcl %in% genes_shared)]))
stopifnot(subcl %in% genes_shared)

print(base::paste("subcl genes used:", 
                  subcl[which(subcl %in% genes_shared)]))


genes_all <- base::unique(genes_subcl$gene)
poss_subcl <- genes_all[which(genes_all %in% genes_shared)]

print(base::paste("subcl genes possible:", 
                  poss_subcl))

#-------------------------------------------------------------------------------
# load and subset annotation from metadata 

anno_subcl <- utils::read.csv(file = anno_subcl_path, 
                              header = TRUE, 
                              sep = ";", 
                              check.names=FALSE, 
                              stringsAsFactors=FALSE, 
                              as.is=TRUE, 
                              colClasses = "character")

print(head(anno_subcl))

anno_subcl <- anno_subcl[anno_subcl$fraction == fraction_curr,]
anno_subcl <- anno_subcl[anno_subcl$cluster == cluster_curr,]

# obtain the nr of subclusters to be separated by mclust 
nr_subcl <- base::max(anno_subcl$subcluster) 
print(cluster_curr)
print(nr_subcl)

#-------------------------------------------------------------------------------
# perform clustering using mclust and the chosen subclustering genes

data <- SingleCellExperiment::logcounts(sce_clust) # normalized multibatchnorm
print(base::paste("final print", subcl))
data <- t(data[rownames(data) %in% subcl,])
head(data)

clust <- mclust::Mclust(data, 
                        G = nr_subcl, 
                        modelNames = c("EII",
                                       "VII",
                                       "EEI",
                                       "VEI",
                                       "EVI",
                                       "VVI",
                                       "EEE", 
                                       "VEE",
                                       "EVE",
                                       "VVE",
                                       "EEV",
                                       "VEV",
                                       "EVV"))
# all multivariate models for data with observations > variables, except VVV
# VVV doesn't work on the cluster, but locally run models all contained one E 

base::summary(clust)

#-------------------------------------------------------------------------------
# add subcluster nr to SCE object (vector between 1 and nr_subcl)
# so that the annotation can be added

sce_clust$subcluster <- base::as.character(clust$classification)
print(base::unique(sce_clust$subcluster))
print(anno_subcl)
 
# add subcluster identity from anno_subcl table
for(i in 1:nr_subcl){
  sce_clust$annotation_subcluster[
    sce_clust$subcluster == i] <- anno_subcl$annotation_subcluster[
      anno_subcl$subcluster == i]
}
print(base::unique(sce_clust$annotation_subcluster))

# remove subcluster number again (is 1 or 2, could lead to confusion later)
sce_clust$subcluster <- NULL

# add info on subclustering genes into rowData
rowData(sce_clust)$subclustering_genes <- vector(length = nrow(sce_clust))
rowData(sce_clust)$subclustering_genes[
  rowData(sce_clust)$Symbol %in% subcl] <- base::paste0(
    base::unique(sce_clust$annotation_subcluster)[1], 
    "_",
    base::unique(sce_clust$annotation_subcluster)[2])
print(table(rowData(sce_clust)$subclustering_genes))

#-------------------------------------------------------------------------------

# save only the sub-clustered portion of the original SCE object
base::saveRDS(sce_clust, snakemake@output[["sce_output"]])

utils::sessionInfo()