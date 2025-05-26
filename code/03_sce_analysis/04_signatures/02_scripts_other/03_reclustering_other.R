
# recluster datasets from other species using different gene sets
# using different resolutions

# determine random number generator for sample
library(parallel)
RNGkind("L'Ecuyer-CMRG") # using this for parallel is necessary

set.seed(37)

#-------------------------------------------------------------------------------

library(Seurat, quietly = TRUE)
library(SeuratObject, quietly = TRUE)

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD 

# params for re-clustering
cut_off_prop <- snakemake@params[["cut_off_prop"]]
nr_cores <- snakemake@params[["nr_cores"]]

#-------------------------------------------------------------------------------

# seurat object from other species

seu <- base::readRDS(snakemake@input[["seu_input"]])

dataset_curr <- snakemake@wildcards[["dataset"]]
print(dataset_curr)

# determine the current fraction to choose the correct gene sets
datasets_other_hsc <- snakemake@params[["datasets_other_hsc"]]
datasets_other_str <- snakemake@params[["datasets_other_str"]]

if(dataset_curr %in% datasets_other_str){
  fraction_curr <- "str"
}else if(dataset_curr %in% datasets_other_hsc){
  fraction_curr <- "hsc"
}
print(fraction_curr)

#-------------------------------------------------------------------------------

seed1 <- (nchar(dataset_curr) + 123)

#-------------------------------------------------------------------------------

# ensembl data frames for each gene set

sign_paths <- snakemake@input[["ensembl_sign"]]
mark_paths <- snakemake@input[["ensembl_mark"]]
mmms_paths <- snakemake@input[["ensembl_mmms"]]

ensembl_sign <- base::readRDS(sign_paths[[
  base::grep(fraction_curr, sign_paths)]])
ensembl_mark <- base::readRDS(mark_paths[[
  base::grep(fraction_curr, mark_paths)]])
ensembl_mmms <- base::readRDS(mmms_paths[[
  base::grep(fraction_curr, mmms_paths)]])

print(head(ensembl_mmms))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# SUBSET

# extract with ensembl column (ID/Symbol/Species) should be used
ensembl_column_use <- seu@misc$ensembl_column_use
print(ensembl_column_use)

# extract gene sets, keep only genes that are also in the Seurat Object
# conserved signature
conserved_signature_IDs <- base::unique(
  ensembl_sign[,which(colnames(ensembl_sign) == ensembl_column_use)])
conserved_signature_IDs <- conserved_signature_IDs[
  conserved_signature_IDs %in% rownames(seu)]

# conserved markers
conserved_marker_IDs <- base::unique(
  ensembl_mark[,which(colnames(ensembl_mark) == ensembl_column_use)])
conserved_marker_IDs <- conserved_marker_IDs[
  conserved_marker_IDs %in% rownames(seu)]

# all BL6 markers
mmusall_marker_IDs <- base::unique(
  ensembl_mmms[,which(colnames(ensembl_mmms) == ensembl_column_use)])
mmusall_marker_IDs <- mmusall_marker_IDs[
  mmusall_marker_IDs %in% rownames(seu)]


# check length and amount of duplication
print("checking ENSMUS ID duplication")
print(base::table(base::duplicated(ensembl_sign$ENSMUS_ID)))

print("checking species ID duplication")
print(base::table(base::duplicated(
  ensembl_sign[,which(colnames(ensembl_sign) == ensembl_column_use)])))

#-------------------------------------------------------------------------------
# random genes
# all genes with expressed in at least cut_off_prop % of cells
# this is just for visualisation

prop_df <- prop_expressed_total_seu(
  seu = seu, 
  geneset = rownames(seu))

# subset by cut-off proportion
prop_df_sub <- prop_df[prop_df$prop_cells >= cut_off_prop,]
print(nrow(prop_df_sub))
stopifnot(nrow(prop_df_sub) > 5000) # so nr of genes is not too small

gene_pool <- prop_df_sub$gene

# get the same nr of random non-0 genes as there are conserved signature genes
# always generate the same random numbers

set.seed(seed1)
random_features <- gene_pool[
  base::sample(1:length(gene_pool), 
               length(conserved_signature_IDs), # same as cons_sign for vis
               replace = F)]
print(seed1)
print(head(random_features))

set.seed(37)

#-------------------------------------------------------------------------------

# subset by gene sets

# set of conserved signature genes
# by using slot = "count" raw counts will be kept in RNA slot 
seu_sign <- BiocGenerics::subset(
  seu, 
  features = conserved_signature_IDs,
  slot = "count")
seu_sign@misc$all_features_reclustering <- conserved_signature_IDs
seu_sign@misc$used_genes <- "conserved_signature"
seu_sign@misc$nr_genes_used <- length(conserved_signature_IDs)

# set of conserved marker genes 
seu_mark <- BiocGenerics::subset(
  seu, 
  features = conserved_marker_IDs,
  slot = "count")
seu_mark@misc$all_features_reclustering <- conserved_marker_IDs
seu_mark@misc$used_genes <- "conserved_markers"
seu_mark@misc$nr_genes_used <- length(conserved_marker_IDs)

# set of all BL6 marker genes 
seu_mmms <- BiocGenerics::subset(
  seu, 
  features = mmusall_marker_IDs,
  slot = "count")
seu_mmms@misc$all_features_reclustering <- mmusall_marker_IDs
seu_mmms@misc$used_genes <- "mmusall_markers"
seu_mmms@misc$nr_genes_used <- length(mmusall_marker_IDs)

# random genes
seu_rand <- BiocGenerics::subset(
  seu, 
  features = random_features,
  slot = "count")
seu_rand@misc$all_features_reclustering <- random_features
seu_rand@misc$used_genes <- "random_features"
seu_rand@misc$nr_genes_used <- length(random_features)

print(seu_sign@misc$nr_genes_used)
print(seu_mark@misc$nr_genes_used)
print(seu_mmms@misc$nr_genes_used)
print(seu_rand@misc$nr_genes_used)

#-------------------------------------------------------------------------------

dim(seu_sign)
dim(seu_mark)
dim(seu_mmms)
dim(seu_rand)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# RECLUSTER

# go through *standard* Seurat pipeline 
# starts from raw counts 

# to choose optimal cluster resolution later
resolution_list <- as.list(seq(0.1, 1.3, by = 0.05))

print("nr_cores")
print(nr_cores)
print("cut_off_prop")
print(cut_off_prop)

# run the standard pipeline for each resolution and each subsetted dataset                           
seu_sign_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, # own function
  seu = seu_sign,
  data_use = seu_sign@misc$data_use,
  calc_umap = TRUE,
  features = conserved_signature_IDs,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE)

# for testing without mclapply
# seu_sign_list_reclustered <- lapply(
#   X = resolution_list,
#   FUN = standard_seu_pipeline, # from source, own function
#   seu = seu_sign,
#   data_use = seu_sign@misc$data_use,
#   calc_umap = TRUE,
#   features = conserved_signature_IDs)

names(seu_sign_list_reclustered) <- base::as.character(unlist(resolution_list))
print(seu_sign_list_reclustered[[1]]@misc)

print("done signature gene reclustering")

seu_mark_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, 
  seu = seu_mark,
  data_use = seu_mark@misc$data_use,
  calc_umap = TRUE,
  features = conserved_marker_IDs,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE)
names(seu_mark_list_reclustered) <- base::as.character(unlist(resolution_list))

print("done conserved marker gene reclustering")

seu_mmms_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, 
  seu = seu_mmms,
  data_use = seu_mmms@misc$data_use,
  calc_umap = TRUE,
  features = mmusall_marker_IDs,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE)
names(seu_mmms_list_reclustered) <- base::as.character(unlist(resolution_list))

print("done BL6 marker gene reclustering")

seu_rand_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, 
  seu = seu_rand,
  data_use = seu_rand@misc$data_use,
  calc_umap = TRUE,
  features = random_features,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE)
names(seu_rand_list_reclustered) <- base::as.character(unlist(resolution_list))

print("done random gene reclustering")

#-------------------------------------------------------------------------------

export_list <- list("seu_sign" = seu_sign_list_reclustered,
                    "seu_mark" = seu_mark_list_reclustered,
                    "seu_mmms" = seu_mmms_list_reclustered,
                    "seu_rand" = seu_rand_list_reclustered)

#-------------------------------------------------------------------------------

base::saveRDS(export_list, snakemake@output[["seu_output"]])

utils::sessionInfo()
