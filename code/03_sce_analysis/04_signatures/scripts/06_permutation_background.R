
# perform permutation tests on reclustering scores as chosen


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

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD OBJECTS

# params for re-clustering and permutation test

cut_off_counts <- snakemake@params[["cut_off_counts"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

# which conservation level to use
cons_level_use <- snakemake@params[["cons_level_use"]]

#-------------------------------------------------------------------------------

# decide which gene sets to use for the current dataset (based on fraction)
dataset_curr <- snakemake@wildcards[["dataset"]]

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

#-------------------------------------------------------------------------------

# load ensembl IDs of specific gene set to be tested from correct fraction
# this can be signature genes or all BL6 marker genes
ensembl_paths <- snakemake@input[["ensembl_paths"]]
print(ensembl_paths)

# this is the correct one for conserved_signature_genes
ensembl_df <- base::readRDS(ensembl_paths[[
  which(base::grepl(fraction_curr, ensembl_paths))]])
print(head(ensembl_df))

#-------------------------------------------------------------------------------

# load the datafram which contains info on which resolution to use for
# re-clustering
resolution_df_path <- snakemake@params[["resolution_df"]] 

resolution_df <- utils::read.csv(file = resolution_df_path, 
                                 header = TRUE, 
                                 sep = ";", 
                                 check.names=FALSE, 
                                 stringsAsFactors=FALSE, 
                                 as.is=TRUE, 
                                 colClasses = "character")

resolution_df <- resolution_df[resolution_df$dataset == dataset_curr,]

# always use the same resolution as for the original gene set
resl <- resolution_df$resolution[
  resolution_df$conservation_level == "random_features"]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE

# extract gene set to be tested, keep only genes that are also in the 
# Seurat Object 
ensembl_column_use <- seu_preprocessed@misc$ensembl_column_use
print(dataset_curr)
print("ensembl_column_use")
print(ensembl_column_use)

test_IDs <- base::unique(
  ensembl_df[,which(colnames(ensembl_df) == ensembl_column_use)])
test_IDs <- test_IDs[
  test_IDs %in% rownames(seu_preprocessed)]

# get the number of conserved signature genes, the same number of random genes 
# will be used for the permutation test
nr_recl_genes <- length(test_IDs)
print("nr_recl_genes")
print(nr_recl_genes)

#-------------------------------------------------------------------------------

# get only genes that have a count of at least n = cut_off_counts
# test genes can be randomly sampled
print(base::summary(rowSums(seu_preprocessed@assays$RNA$counts)))

gene_pool <- rownames(seu_preprocessed)[
  which(rowSums(seu_preprocessed@assays$RNA$counts) > cut_off_counts)]
print(length(gene_pool))

# subset seurat object to these genes as the other genes won't be used
seu_preprocessed_sub <- BiocGenerics::subset(seu_preprocessed,
                                             features = gene_pool,
                                             slot = "count")
print(dim(seu_preprocessed_sub))
print(base::summary(rowSums(seu_preprocessed_sub@assays$RNA$counts)))

#-------------------------------------------------------------------------------

# generate i=iterations random sets of the same length as there are test genes
# always generate the same random numbers

iteration_df <- base::data.frame(row.names = c(1:nr_recl_genes))

set.seed(37)
base::RNGkind("L'Ecuyer-CMRG")
for(i in 1:iterations){
  iteration_df[,i] <- base::sample(1:length(gene_pool), 
                                   nr_recl_genes, 
                                   replace = FALSE)
  colnames(iteration_df)[i] <- i
}
print(iteration_df[1:10,1:2])

#-------------------------------------------------------------------------------
# run standard seurat re-clustering pipeline (own function) i times

print("nr_cores")
print(nr_cores)
print("cut_off_counts")
print(cut_off_counts)
print("resl")
print(resl)
print("starting iteration")

res_df_list <- parallel::mclapply(X = as.list(c(1:iterations)),
                                  FUN = random_reclustering_scores,
                                  seu = seu_preprocessed_sub,
                                  assay_use = "RNA",
                                  iteration_df = iteration_df,
                                  resolution = resl,
                                  mc.preschedule = TRUE,
                                  mc.cores = nr_cores,
                                  mc.silent = TRUE)

# res_df_list <- lapply(X = as.list(c(1:iterations)),
#                                   FUN = random_reclustering_scores,
#                                   seu = seu_preprocessed_sub,
#                                   assay_use = "RNA",
#                                   iteration_df = iteration_df,
#                                   resolution = resl)

score_df <- dplyr::bind_rows(res_df_list)
print(head(score_df))

#-------------------------------------------------------------------------------

base::saveRDS(score_df, snakemake@output[["perm_score_df"]])

utils::sessionInfo()
