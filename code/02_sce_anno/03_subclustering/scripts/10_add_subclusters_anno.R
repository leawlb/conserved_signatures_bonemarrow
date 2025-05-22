#-------------------------------------------------------------------------------
# add annotation to subclustered 
# recalculate nicer looking UMAP after finalizing cell type annotation

RNGkind("L'Ecuyer-CMRG") 
set.seed(37)

library(scran, quietly = TRUE)
library(scater, quietly = TRUE)

#-------------------------------------------------------------------------------
# load objects

# original SCE object with only cluster annotations
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

fraction_curr <- snakemake@wildcards[["fraction"]]

sce_subcl_path <- snakemake@input[["sce_subcl"]]
anno_final_path <- snakemake@input[["anno_final"]]

nr_hvgs <- snakemake@params[["nr_hvgs"]]
seeds_umap <- snakemake@params[["seeds_umap"]]

seed <- seeds_umap[[fraction_curr]]

#-------------------------------------------------------------------------------
# load subclustered SCEs
# these are ALL subclustered SCE subsets, from both fractions

sce_subcl_list <- list()
for(i in 1:length(sce_subcl_path)){
  sce_subcl_list[[i]] <- base::readRDS(file = sce_subcl_path[[i]])
}

# subset to current fraction
sce_subcl_list <- lapply(sce_subcl_list, function(sce){
  if(fraction_curr %in% sce$Fraction_ID){
    return(sce)
  }
})

# remove any empty list items
sce_subcl_list <- sce_subcl_list[lengths(sce_subcl_list) != 0]
print(sce_subcl_list)

#-------------------------------------------------------------------------------
# load final annotation for clusters and subclusters = cell types

final_anno <- utils::read.csv(file = anno_final_path, 
                              header = TRUE, 
                              sep = ";", 
                              check.names=FALSE, 
                              stringsAsFactors=FALSE, 
                              as.is=TRUE, 
                              colClasses = "character")
final_anno <- final_anno[final_anno$fraction == fraction_curr,]

print(head(final_anno)) 

#-------------------------------------------------------------------------------
# add info

sce$annotation_subcluster <- vector(length = ncol(sce))
rowData(sce)$subclustering_genes <- vector(length = nrow(sce))

# for each subclustered SCE object
for(i in 1:length(sce_subcl_list)){
  sce_temp <- sce_subcl_list[[i]]

  # check that colnames = barcodes match 
  print(base::table(is.na(base::match(colnames(sce_temp), colnames(sce)))))
  print(base::table(is.na(base::match(colnames(sce), colnames(sce_temp)))))
  
  # add subcluster annotation from the subclustered to the original sce
  # this includes ONLY the subcluster ANNOTATION = cell type name
  sce$annotation_subcluster[
    base::match(colnames(sce_temp), 
                colnames(sce))] <- sce_temp$annotation_subcluster
  
  # get info on which genes were used for subclustering which cell types
  subcl <-  rowData(sce_temp)$Symbol[
    rowData(sce_temp)$subclustering_genes != FALSE] 
  
  # add info to sce
  rowData(sce)$subclustering_genes[
    rowData(sce)$Symbol %in% subcl] <- rowData(sce_temp)$subclustering_genes[
      rowData(sce_temp)$subclustering_genes != FALSE][1]
}

# fill empty slots with cluster names (cluster annotation = cell type names)
# these are cell types that were not sub-clustered
sce$annotation_subcluster[
  sce$annotation_subcluster == FALSE] <- unfactor(sce$annotation_cluster[
    sce$annotation_subcluster == FALSE])

print(base::unique(sce$annotation_subcluster))

#-------------------------------------------------------------------------------
# add final annotation into "celltypes" + category
# also add cluster_after_sub

sce$celltypes <- sce$annotation_subcluster

sce$category <- vector(length = ncol(sce))
sce$cluster_after_sub <- vector(length = ncol(sce))
for(i in 1:length(base::unique(sce$celltypes))){
  
  ct <- base::unique(sce$celltypes)[i]
  
  sce$category[sce$celltypes == ct] <- final_anno$category[
    final_anno$celltypes == ct]
  
  sce$cluster_after_sub[sce$celltypes == ct] <- i
}
sce$category <- factor(
  sce$category, levels = base::unique(final_anno$category))

sce$cluster_after_sub <- factor(
  sce$cluster_after_sub, 
  levels = c(1:max(base::as.numeric(sce$cluster_after_sub))))

#-------------------------------------------------------------------------------
# test

print("TEST")

print(base::sort(base::unique(final_anno$celltypes)))
print(base::sort(base::unique(sce$celltypes)))

stopifnot(!is.na(base::match(base::unique(final_anno$celltypes),
                             base::unique(sce$celltypes))))
stopifnot(!is.na(base::match(base::unique(sce$celltypes), 
                             base::unique(final_anno$celltypes))))

print(final_anno$celltypes)

print(base::table(sce$celltypes, sce$annotation_subcluster))
print(base::table(sce$celltypes, sce$annotation_cluster))

sce$celltypes <- factor(sce$celltypes, levels = final_anno$celltypes)

print(levels(sce$category))
print(levels(sce$celltypes))

stopifnot(!is.na(sce$category))
stopifnot(!is.na(sce$celltypes))
stopifnot(!is.na(sce$cluster_after_sub))
stopifnot(!is.na(sce$annotation_subcluster))

print(base::table(rowData(sce)$subclustering_genes))

#-------------------------------------------------------------------------------
# re-calculate UMAP coordinates for better looking plots
# now that some clusters have been removed and all cell types are annotated

set.seed(seed)
print(seed)

gene_var <- scran::modelGeneVar(sce)
hvgs <- scran::getTopHVGs(gene_var, n = nr_hvgs)

# use batch-corrected PC coordinates in "PCA"
sce <- scater::runUMAP(sce, dimred = "PCA", subset_row = hvgs)
colnames(SingleCellExperiment::reducedDims(sce)$UMAP) <- c("X1", "X2")

#-------------------------------------------------------------------------------

base::saveRDS(sce, snakemake@output[["sce_output"]])

base::saveRDS(sce, snakemake@output[["sce_output_pretty"]])

utils::sessionInfo()