
# EXPLANATION TODO


#-------------------------------------------------------------------------------
set.seed(37)
# base::RNGkind("L'Ecuyer-CMRG") is downstream

library(Seurat, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(BiocGenerics, quietly = TRUE)
library(S4Vectors, quietly = TRUE)

# for using mclapply
library(parallel, quietly = TRUE)

# for calculating scores
library(mclust, quietly = TRUE)
library(mcclust, quietly = TRUE)
library(bluster, quietly = TRUE)
library(dendextend, quietly = TRUE)

source(snakemake@params[["reclustering_functions"]])
#source("repositories/Interspecies_BM_phd/code/source/sce_functions_reclustering.R")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD OBJECTS

# params for re-clustering and permutation test

# cut_off_counts <- 5
# nr_cores <- 1
# iterations <- 2

cut_off_counts <- snakemake@params[["cut_off_counts"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

#-------------------------------------------------------------------------------

# decide which gene sets to use for the current dataset (based on fraction)
dataset_curr <- snakemake@wildcards[["dataset"]]

# dataset_curr <- "li_all_stromal"

if(dataset_curr %in% c("ts_all_stromal", "li_all_stromal")){
  fraction_curr <- "str"
}else if(dataset_curr %in% c("ts_bone_marrow", 
                             "ts_hscs_progenitors",
                             "nmr_sorted_hspc",
                             "zeb_all_hspc")){
  fraction_curr <- "hsc"
}

#-------------------------------------------------------------------------------

# load the seurat object/dataset to be tested
seu_preprocessed <- base::readRDS(snakemake@input[["seu_preprocessed"]])

#seu_preprocessed <- base::readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/prepared/li_all_stromal")

#-------------------------------------------------------------------------------

# load ensembl IDs of specific gene set to be tested from correct fraction
# this can be signature genes 
ensembl_sign <- snakemake@input[["ensembl_sign"]]
ensembl_mark <- snakemake@input[["ensembl_mark"]]
ensembl_mmms <- snakemake@input[["ensembl_mmms"]]

#ensembl_sign <- c("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/02_endf/ensembl_sign_str")
#ensembl_mark <- c("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/02_endf/ensembl_mark_str")
#ensembl_mmms <- c("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/02_endf/ensembl_mmms_str")

ensembl_sign_df <- base::readRDS(ensembl_sign[[
  which(base::grepl(fraction_curr, ensembl_sign))]])
ensembl_mark_df <- base::readRDS(ensembl_mark[[
  which(base::grepl(fraction_curr, ensembl_mark))]])
ensembl_mmms_df <- base::readRDS(ensembl_mmms[[
  which(base::grepl(fraction_curr, ensembl_mmms))]])

#-------------------------------------------------------------------------------

# load the dataframe which contains info on which resolution to use for
# re-clustering
resolution_df_path <- snakemake@params[["resolution_df"]] 
#resolution_df_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/reclustering_other_resolution.txt"

resolution_df <- utils::read.csv(file = resolution_df_path, 
                                 header = TRUE, 
                                 sep = ";", 
                                 check.names=FALSE, 
                                 stringsAsFactors=FALSE, 
                                 as.is=TRUE, 
                                 colClasses = "character")

resolution_df <- resolution_df[resolution_df$dataset == dataset_curr,]

# get resolution for each gene set to be permuted
resl_mark <- resolution_df$resolution[
  resolution_df$conservation_level == "conserved_markers"]

resl_mmms <- resolution_df$resolution[
  resolution_df$conservation_level == "mmusall_markers"]

print(resl_mark)
print(resl_mmms)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE

# extract conserved signature gene set
ensembl_column_use <- seu_preprocessed@misc$ensembl_column_use
print(dataset_curr)
print("ensembl_column_use")
print(ensembl_column_use)

sign_IDs <- base::unique(
  ensembl_sign_df[,which(colnames(ensembl_sign_df) == ensembl_column_use)])
sign_IDs <- sign_IDs[sign_IDs %in% rownames(seu_preprocessed)]

# conserved markers
mark_IDs <- base::unique(
  ensembl_mark_df[,which(colnames(ensembl_mark_df) == ensembl_column_use)])
mark_IDs <- mark_IDs[mark_IDs %in% rownames(seu_preprocessed)]

# all BL6 markers
mmms_IDs <- base::unique(
  ensembl_mmms_df[,which(colnames(ensembl_mmms_df) == ensembl_column_use)])
mmms_IDs <- mmms_IDs[mmms_IDs %in% rownames(seu_preprocessed)]

# get the number of random genes that need to be added to 
# conserved_signature_IDs for each gene set to be permuted
nr_sign <- length(sign_IDs) 
nr_random_mark <- length(mark_IDs) - nr_sign
nr_random_mmms <- length(mmms_IDs) - nr_sign

print(nr_sign)
print(nr_random_mark)
print(length(mark_IDs))
print(nr_random_mmms)
print(length(mmms_IDs))

#-------------------------------------------------------------------------------

# get position of conserved signature genes 
# POSITIONS FROM ORIGINAL SEU_PREPROCESSED

SIGN_POSITIONS <- which(rownames(seu_preprocessed) %in% sign_IDs)
print(length(SIGN_POSITIONS))

# separate seurat object into two, one with only signature IDs, one without
# the one without will be used as pool for getting random genes
seu_sign <- BiocGenerics::subset(seu_preprocessed,
                                 features = sign_IDs,
                                 slot = "count")

non_seu_sign <- rownames(seu_preprocessed)[-SIGN_POSITIONS]
seu_pool <- BiocGenerics::subset(seu_preprocessed,
                                 features = non_seu_sign,
                                 slot = "count")

print(nrow(seu_preprocessed))
print(nrow(seu_sign))
print(nrow(seu_pool))

#-------------------------------------------------------------------------------

# get only genes that have a count of at least n = cut_off_counts
print(base::summary(rowSums(seu_pool@assays$RNA$counts)))

gene_pool <- rownames(seu_pool)[
  which(rowSums(seu_pool@assays$RNA$counts) > cut_off_counts)]
print(length(gene_pool))

# subset seurat object to these genes as the other genes won't be used
seu_pool <- BiocGenerics::subset(seu_pool,
                                 features = gene_pool,
                                 slot = "count")

# the positions of genes in SEU that are allowed as pool for random genes
# POSITIONS FROM ORIGINAL SEU_PREPROCESSED
POOL_POSITIONS <- which(rownames(seu_preprocessed) %in% rownames(seu_pool))

print(length(POOL_POSITIONS))

# this is the pool of genes from which random genes can be chosen
# e.g. gene that are not part of conserved signatures genes
# and that are expressed > cut_off_counts

#-------------------------------------------------------------------------------

# generate i = iterations random sets of genes at the required number
# from the pool of random genes (original positions in SEU)
# always generate the same random numbers
set.seed(37)
base::RNGkind("L'Ecuyer-CMRG")

iteration_df_mark <- base::data.frame(row.names = c(1:nr_random_mark))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED subsetted to a smaller list 
  # basically, subset the vector of pool positions in random positions 
  iteration_df_mark[,i] <- POOL_POSITIONS[base::sample(1:length(gene_pool), 
                                                       nr_random_mark, 
                                                       replace = FALSE)]
  stopifnot(iteration_df_mark[,i] %in% POOL_POSITIONS)
  stopifnot(!iteration_df_mark[,i] %in% SIGN_POSITIONS)
  
  colnames(iteration_df_mark)[i] <- i
}
print(iteration_df_mark[1:10,1:2])
print(dim(iteration_df_mark))

#-------------------------------------------------------------------------------
set.seed(37)
base::RNGkind("L'Ecuyer-CMRG")

iteration_df_mmms <- base::data.frame(row.names = c(1:nr_random_mmms))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED subsetted to a smaller list 
  iteration_df_mmms[,i] <- POOL_POSITIONS[base::sample(1:length(gene_pool), 
                                                       nr_random_mmms, 
                                                       replace = FALSE)]
  stopifnot(iteration_df_mmms[,i] %in% POOL_POSITIONS)
  stopifnot(!iteration_df_mmms[,i] %in% SIGN_POSITIONS)
  
  colnames(iteration_df_mmms)[i] <- i
}
print(iteration_df_mmms[1:10,1:2])
print(dim(iteration_df_mmms))

#-------------------------------------------------------------------------------
# generate a fitting data frame to add signature gene positions to the 
# random positions from the allowed pool

iteration_df_sign <- base::data.frame(row.names = c(1:nr_sign))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED
  iteration_df_sign[,i] <- SIGN_POSITIONS
  colnames(iteration_df_sign)[i] <- i
}
print(iteration_df_sign[1:10,1:2])
print(dim(iteration_df_sign))

#-------------------------------------------------------------------------------

iteration_df_mark <- base::rbind(iteration_df_sign, iteration_df_mark)
iteration_df_mmms <- base::rbind(iteration_df_sign, iteration_df_mmms)

# make sure there are no duplicated positions as failsave
# ALL POSITIONS FROM ORIGINAL SEU_PREPROCESSED
for(i in 1:iterations){
  stopifnot(!duplicated(iteration_df_mark[,i]))
  stopifnot(!duplicated(iteration_df_mmms[,i]))
}

#-------------------------------------------------------------------------------
# run standard seurat re-clustering pipeline (own function) i times

print("nr_cores")
print(nr_cores)
print("cut_off_counts")
print(cut_off_counts)
print("resl")
print(resl_mark)
print(resl_mmms)
print("starting iteration mark")

res_df_list_mark <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = random_reclustering_scores,
  seu = seu_preprocessed,
  assay_use = "RNA",
  iteration_df = iteration_df_mark,
  resolution = resl_mark,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

# res_df_list_mark <- lapply(
#   X = as.list(c(1:iterations)),
#   FUN = random_reclustering_scores,
#   seu = seu_preprocessed,
#   assay_use = "RNA",
#   iteration_df = iteration_df_mark,
#   resolution = resl_mark)

score_df_mark <- dplyr::bind_rows(res_df_list_mark)
print(head(score_df_mark))

print("starting iteration mmms")

res_df_list_mmms <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = random_reclustering_scores,
  seu = seu_preprocessed,
  assay_use = "RNA",
  iteration_df = iteration_df_mmms,
  resolution = resl_mmms,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

score_df_mmms <- dplyr::bind_rows(res_df_list_mmms)
print(head(score_df_mmms))

#-------------------------------------------------------------------------------

base::saveRDS(score_df_mark, snakemake@output[["perm_score_df_mark"]])
base::saveRDS(score_df_mmms, snakemake@output[["perm_score_df_mmms"]])

utils::sessionInfo()
