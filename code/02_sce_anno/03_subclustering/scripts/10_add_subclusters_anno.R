#-------------------------------------------------------------------------------

library(scran, quietly = TRUE)
library(scater, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
sce_subcl_path <- snakemake@input[["sce_subcl"]]
fraction_curr <- snakemake@wildcards[["fraction"]]

nr_hvgs <- snakemake@params[["nr_hvgs"]]
seeds_umap <- snakemake@params[["seeds_umap"]]

if(fraction_curr == "hsc"){
  seed <- seeds_umap[["hsc"]]
}else if(fraction_curr == "str"){
  seed <- seeds_umap[["str"]]
}

#-------------------------------------------------------------------------------

print(sce)

sce_subcl_list <- list()
for(i in 1:length(sce_subcl_path)){
  sce_subcl_list[[i]] <- readRDS(file = sce_subcl_path[[i]])
}

sce_subcl_list <- lapply(sce_subcl_list, function(sce){
  if(fraction_curr %in% sce$Fraction_ID){
    return(sce)
  }
})

sce_subcl_list <- sce_subcl_list[lengths(sce_subcl_list) != 0]
print(sce_subcl_list)

#-------------------------------------------------------------------------------
# load final annotation for clusters and subclusters = cell types

sce$annotation_subcluster <- vector(length = ncol(sce))
sce$subcluster <- vector(length = ncol(sce))

final_anno <- read.csv(file = snakemake@input[["anno_final"]], 
                       header = TRUE, 
                       sep = ";", 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       as.is=TRUE, 
                       colClasses = "character")
final_anno <- final_anno[final_anno$fraction == fraction_curr,]

head(final_anno)

#-------------------------------------------------------------------------------
# add more info

for(i in 1:length(sce_subcl_list)){
  sce_temp <- sce_subcl_list[[i]]

  print(table(is.na(match(colnames(sce_temp), colnames(sce)))))
  print(table(is.na(match(colnames(sce), colnames(sce_temp)))))
  
  sce$annotation_subcluster[
    match(colnames(sce_temp), colnames(sce))] <- sce_temp$annotation_subcluster
  sce$subcluster[
    match(colnames(sce_temp), colnames(sce))] <- sce_temp$subcluster
}

# fill empty slots with cluster names
sce$annotation_subcluster[
  sce$annotation_subcluster == FALSE] <- unfactor(sce$annotation_cluster[
    sce$annotation_subcluster == FALSE])
sce$subcluster[
  sce$subcluster == FALSE] <- sce$cluster_louvain[sce$subcluster == FALSE]

print(unique(sce$annotation_subcluster))
print(unique(sce$subcluster))

#-------------------------------------------------------------------------------
# add final annotation into "celltypes" 

sce$celltypes <- sce$annotation_subcluster

print(sort(unique(final_anno$celltypes)))
print(sort(unique(sce$celltypes)))

stopifnot(!is.na(match(unique(final_anno$celltypes), unique(sce$celltypes))))
stopifnot(!is.na(match(unique(sce$celltypes), unique(final_anno$celltypes))))

print(final_anno$celltypes)
sce$celltypes <- factor(sce$celltypes, levels = final_anno$celltypes)

#-------------------------------------------------------------------------------
# add category 

sce$category <- vector(length = ncol(sce))

for(i in unique(sce$celltypes)){
  sce$category[sce$celltypes == i] <- final_anno$category[
    final_anno$celltypes == i]
}
sce$category <- factor(sce$category, levels = unique(final_anno$category))

print(levels(sce$category))
print(levels(sce$celltypes))

stopifnot(!is.na(sce$category))
stopifnot(!is.na(sce$celltypes))

#-------------------------------------------------------------------------------
# re-calculate UMAP coordinates for better looking plots
# now that some cell types are missing

gene_var <- scran::modelGeneVar(sce)
hvgs_for_batch_correction <- scran::getTopHVGs(gene_var, n = nr_hvgs)
hvgs <- scran::getTopHVGs(gene_var, n = nr_hvgs)

set.seed(seed)

# use batch-corrected PC coordinates in "PCA
sce <- scater::runUMAP(sce, dimred = "PCA", subset_row = hvgs)
colnames(reducedDims(sce)$UMAP) <- c("X1", "X2")

set.seed(37)

#-------------------------------------------------------------------------------

saveRDS(sce, snakemake@output[["sce_output"]])

sessionInfo()