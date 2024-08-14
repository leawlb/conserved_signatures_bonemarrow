
# recluster datasets from other species using different gene sets

library(Seurat, quietly = TRUE)
library(SeuratObject, quietly = TRUE)
library(parallel, quietly = TRUE)

set.seed(37)
base::RNGkind("L'Ecuyer-CMRG") # for mclapply, generation of same random numbers

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD 

# params for re-clustering
cut_off_counts <- snakemake@params[["cut_off_counts"]]
nr_cores <- snakemake@params[["nr_cores"]]

#-------------------------------------------------------------------------------

# seurat object from other species

seu <- base::readRDS(snakemake@input[["seu_input"]])

dat_curr <- snakemake@wildcards[["dataset"]]
print(dat_curr)

# determine the current fraction to choose the correct gene sets
if(dat_curr %in%  c("ts_all_stromal", "li_all_stromal")){
  fraction_curr <- "str"
}else if(dat_curr %in%  c("ts_bone_marrow", 
                          "ts_hscs_progenitors",
                          "nmr_sorted_hspc",
                          "zeb_all_hspc")){
  fraction_curr <- "hsc"
}
print(fraction_curr)

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
print("checking ENSMUS ID")
print(base::table(base::duplicated(ensembl_sign$ENSMUS_ID)))

print("checking species ID")
print(base::table(base::duplicated(
  ensembl_sign[,which(colnames(ensembl_sign) == ensembl_column_use)])))

# random genes
# all genes with non-zero expression (above a defined threshold)
non_zero_features <- rownames(seu@assays$RNA)[
  BiocGenerics::rowSums(seu@assays$RNA) > cut_off_counts]

# get the same nr of random non-0 genes as there are conserved signature genes
random_features <- non_zero_features[
  base::sample(1:length(non_zero_features), 
               length(conserved_signature_IDs),
               replace = F)]

#-------------------------------------------------------------------------------

# subset by gene sets

# set of conserved signature genes
# by using slot = "count" raw counts will be kept in RNA slot 
seu_sign <- BiocGenerics::subset(
  seu, 
  features = conserved_signature_IDs,
  slot = "count")
seu_sign@misc$all_features_subclustering <- conserved_signature_IDs
seu_sign@misc$used_genes <- "conserved_signature"

# set of conserved marker genes 
seu_mark <- BiocGenerics::subset(
  seu, 
  features = conserved_marker_IDs,
  slot = "count")
seu_mark@misc$all_features_subclustering <- conserved_marker_IDs
seu_mark@misc$used_genes <- "conserved_markers"

# set of all BL6 marker genes 
seu_mmms <- BiocGenerics::subset(
  seu, 
  features = mmusall_marker_IDs,
  slot = "count")
seu_mmms@misc$all_features_subclustering <- mmusall_marker_IDs
seu_mmms@misc$used_genes <- "mmusall_markers"

# random genes
seu_rand <- BiocGenerics::subset(
  seu, 
  features = random_features,
  slot = "count")
seu_rand@misc$all_features_subclustering <- random_features
seu_rand@misc$used_genes <- "random_features"

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
#resolution_list <- as.list(c(0.5, 0.1))

print(nr_cores)
print(cut_off_counts)
print(resolution_list)

#nr_cores <- 1
# run the standard pipeline for each resolution and each subsetted dataset                           
seu_sign_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, # from source, own function
  seu = seu_sign,
  assay_use = "RNA",
  calc_umap = TRUE,
  features = conserved_signature_IDs,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE)
names(seu_sign_list_reclustered) <- base::as.character(unlist(resolution_list))

print("done signature gene reclustering")

seu_mark_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, 
  seu = seu_mark,
  assay_use = "RNA",
  calc_umap = TRUE,
  features = conserved_marker_IDs,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE)
names(seu_mark_list_reclustered) <- base::as.character(unlist(resolution_list))

print("done conserved marker gene reclustering")

seu_mmms_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, 
  seu = seu_mmms,
  assay_use = "RNA",
  calc_umap = TRUE,
  features = mmusall_marker_IDs,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE)
names(seu_mmms_list_reclustered) <- base::as.character(unlist(resolution_list))

print("done BL6 marker gene reclustering")

seu_rand_list_reclustered <- mclapply(
  X = resolution_list,
  FUN = standard_seu_pipeline, 
  seu = seu_rand,
  assay_use = "RNA",
  calc_umap = TRUE,
  features = random_features,
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE)
names(seu_rand_list_reclustered) <- base::as.character(unlist(resolution_list))

print("done random gene reclustering")

#-------------------------------------------------------------------------------

export_list <- list("seu_sign" = seu_sign_list_reclustered,
                    "seu_mark" = seu_mark_list_reclustered,
                    "seu_mmms" = seu_mmms_list_reclustered,
                    "seu_rand" = seu_rand_list_reclustered)

base::saveRDS(export_list, snakemake@output[["seu_output"]])

utils::sessionInfo()
