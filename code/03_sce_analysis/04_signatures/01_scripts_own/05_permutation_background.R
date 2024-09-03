
#-------------------------------------------------------------------------------

# permutation of own dataset reclustering

set.seed(37)
base::RNGkind("L'Ecuyer-CMRG")

#-------------------------------------------------------------------------------

library(scater, quietly = TRUE)
library(scran, quietly = TRUE)
library(bluster, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(parallel, quietly = TRUE)
library(mclust, quietly = TRUE)
library(mcclust, quietly = TRUE)
library(dendextend, quietly = TRUE)

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# load objects

#-------------------------------------------------------------------------------

# dataset 
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

#-------------------------------------------------------------------------------
# params for re-clustering

k_graph_list <- snakemake@params[["k_graph_list"]]
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]
resolution_louvain <- resolution_louvain_list[[fraction_curr]]

print(k_graph)
print(resolution_louvain)

#-------------------------------------------------------------------------------
# cell types that were not used to extract signature are excluded because
# they cannot be separated after excluding them
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)

#-------------------------------------------------------------------------------
# params for permutation test

cut_off_counts <- snakemake@params[["cut_off_counts"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

#-------------------------------------------------------------------------------
# which conservation level to use
cons_level_use <- snakemake@params[["cons_level_use"]]
print(cons_level_use)

#-------------------------------------------------------------------------------
# ensembl IDs of specific gene set to be tested from correct fraction
ensembl <- snakemake@input[["ensembl"]]

ensembl_df <- base::readRDS(ensembl)
print(head(ensembl_df))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE

#-------------------------------------------------------------------------------
# extract gene set to be tested
ensembl_column_use <- "MMUS_SYMBOL"
print("ensembl_column_use")
print(ensembl_column_use)

test_IDs <- base::unique(
  ensembl_df[,which(colnames(ensembl_df) == ensembl_column_use)])

# get the number of genes to be tested
# the same number of random genes will be used for the permutation test
nr_recl_genes <- length(test_IDs)
print("nr_recl_genes")
print(nr_recl_genes)

#-------------------------------------------------------------------------------

# get only genes that have a count of at least n = cut_off_counts
# test genes can be randomly sampled
# batch normalised
print(base::summary(rowSums(SummarizedExperiment::assays(sce)$logcounts)))

gene_pool <- rownames(sce)[
  which(rowSums(SummarizedExperiment::assays(sce)$logcounts) > cut_off_counts)]
print(length(gene_pool))

# subset object to these genes as the other genes won't be used
# subset GENES
sce_sub <- sce[which(rownames(sce) %in% gene_pool),]

print(dim(sce_sub))
print(nrow(sce_sub))
print(base::summary(rowSums(SummarizedExperiment::assays(sce_sub)$logcounts)))

#-------------------------------------------------------------------------------

# generate i=iterations random sets of the same length as there are test genes
# always generate the same random numbers
# these are the random positions of genes that will be used for re-clustering

iteration_df <- base::data.frame(row.names = c(1:nr_recl_genes))

for(i in 1:iterations){
  iteration_df[,i] <- base::sample(1:length(gene_pool), 
                                   nr_recl_genes, 
                                   replace = FALSE)
  colnames(iteration_df)[i] <- i
}
print(iteration_df[1:10,1:2])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ITERATIONS

print("nr_cores")
print(nr_cores)
print("iterations")
print(iterations)
print("cut_off_counts")
print(cut_off_counts)
print("k_graph")
print(k_graph)
print("resolution_louvain")
print(resolution_louvain)

# run re-clustering and score function at once for i iterations
res_df_list <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = permuting_reclustering_scores_sce,
  sce = sce_sub,
  iteration_df = iteration_df,
  resolution_louvain = resolution_louvain,
  k_graph = k_graph,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

# for testing without mclapply
# res_df_list <- lapply(
#   X = as.list(c(1:iterations)),
#   FUN = random_reclustering_scores_orig,
#   sce = sce_sub,
#   iteration_df = iteration_df,
#   k_graph = k_graph,
#   resolution_louvain = resolution_louvain)

score_df <- dplyr::bind_rows(res_df_list)
print(head(score_df))

#-------------------------------------------------------------------------------
base::saveRDS(score_df, snakemake@output[["perm_score_df"]])

utils::sessionInfo()
