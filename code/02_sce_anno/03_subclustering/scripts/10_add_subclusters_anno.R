#-------------------------------------------------------------------------------

library(SingleCellExperiment)

sce <- readRDS(file = snakemake@input[["sce_input"]])
sce_subcl_path <- snakemake@input[["sce_subcl"]]

sce_subcl_list <- list()
for(i in 1:length(sce_subcl_path)){
  sce_subcl_list[[i]] <- readRDS(file = sce_subcl_path[[i]])
}

print(sce_subcl_list)

sce$annotation_subcluster <- vector(length = ncol(sce))
sce$subcluster <- vector(length = ncol(sce))

fraction_curr <- snakemake@wildcards[["fraction"]]

final_anno <- read.csv(file = snakemake@input[["final_annotation"]], 
                       header = TRUE, 
                       sep = ";", 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       as.is=TRUE, 
                       colClasses = "character")
final_anno <- final_anno[final_anno$fraction == fraction_curr,]

head(final_anno)

#-------------------------------------------------------------------------------
# add info
for(i in 1:length(sce_subcl_list)){
  sce_temp <- sce_subcl_list[[i]]
  
  print(sce_temp)
  print(table(is.na(match(colnames(sce_temp), colnames(sce)))))
  print(table(is.na(match(colnames(sce), colnames(sce_temp)))))
  
  unique(sce_temp$annotation_subcluster)
  unique(sce_temp$subcluster)
  
  sce$annotation_subcluster[
    match(colnames(sce_temp), colnames(sce))] <- sce_temp$annotation_subcluster
  sce$subcluster[
    match(colnames(sce_temp), colnames(sce))] <- sce_temp$subcluster
}

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

stopifnot(!is.na(match(final_anno$celltypes, sce$celltypes)))
stopifnot(!is.na(match(sce$celltypes, final_anno$celltypes)))

sce$celltypes <- factor(sce$celltypes, levels = final_anno$celltypes)
print(table(is.na(sce$celltypes)))
#-------------------------------------------------------------------------------
# add category

sce$category <- vector(length = ncol(sce))

for(i in unique(sce$celltypes)){
  print(i)
  sce$category[sce$celltypes == i] <- final_anno$category[final_anno$celltypes == i]
}
sce$category <- factor(sce$category, levels = unique(final_anno$category))

print(sce$category)
print(sce$celltypes)

saveRDS(sce, snakemake@output[["sce_output"]])