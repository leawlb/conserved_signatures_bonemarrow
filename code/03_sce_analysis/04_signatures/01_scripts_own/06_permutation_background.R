
#-------------------------------------------------------------------------------

# permutation of own dataset reclustering

# determine random number generator for sample
# Mersenne-Twister" is default
RNGkind("Mersenne-Twister") 
set.seed(37)

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

source(snakemake@params[["functions_reclustering"]])

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

#-------------------------------------------------------------------------------
# cell types that were not used to extract signature are excluded because
# they cannot be separated after excluding them
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)
print(base::table(sce$celltypes))

sce <- sce[,which(!sce$celltypes %in% cts_exclude)]
print(dim(sce))
print(base::table(sce$celltypes))

#-------------------------------------------------------------------------------
# params for permutation test

cut_off_prop <- snakemake@params[["cut_off_prop"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

#-------------------------------------------------------------------------------
# which conservation level to use
cons_level_use <- snakemake@params[["cons_level_use"]]
print("cons_level_use")
print(cons_level_use)

#-------------------------------------------------------------------------------
# load geneset list
geneset_list <- base::readRDS(snakemake@input[["geneset_list"]])

#-------------------------------------------------------------------------------
# set seeds for sample()

if(fraction_curr == "hsc"){
  seed1 <- 4000
}else if(fraction_curr == "str"){
  seed1 <- 1234
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE

#-------------------------------------------------------------------------------
# extract subclustering genes for control
subclustering_gene_list <- lapply(geneset_list, function(geneset){
  subclustering_genes <- geneset$genes_subclustering
  return(subclustering_genes)
})
subclustering_genes <- base::unique(unlist(subclustering_gene_list))
stopifnot(!subclustering_genes %in% rownames(sce))

# extract gene set to be tested (based on "cons_level_use")
test_IDs_list <- lapply(geneset_list, function(geneset){
  test_ids <- geneset[[which(names(geneset) == cons_level_use)]]
  return(test_ids)
})
test_IDs <- base::unique(unlist(test_IDs_list))
test_IDs <- test_IDs[test_IDs %in% rownames(sce)]

# get the number of genes to be tested
# the same number of random genes will be used for the permutation test
nr_recl_genes <- length(test_IDs)
print("nr_recl_genes")
print(nr_recl_genes)


#-------------------------------------------------------------------------------

# get only genes that are expressed in >= cut_off_prop% of cells

gene_pool <- rownames(sce)
print(length(gene_pool))

# using own function
prop_df <- prop_expressed_total_sce(
  sce = sce, 
  geneset = gene_pool)

# subset by cut-off proportion
prop_df_sub <- prop_df[prop_df$prop_cells >= cut_off_prop,]
print(nrow(prop_df_sub))

# subset gene pool
gene_pool <- gene_pool[which(gene_pool %in% prop_df_sub$gene)]
print(length(gene_pool))

# subset object to these genes as the other genes won't be used
# subset GENES
sce_sub <- sce[which(rownames(sce) %in% gene_pool),]

print(dim(sce))
print(dim(sce_sub))

#-------------------------------------------------------------------------------

# generate n=iterations random sets of the same length as there are test genes
# always generate the same random numbers (RNG)
# these are the random positions of genes that will be used for re-clustering
# - no subclustering genes
# - no genes expressed below cut_off_prop
# - CAN randomly contain genes from gene set to be tested

iteration_df <- base::data.frame(row.names = c(1:nr_recl_genes))

set.seed(seed1)
for(i in 1:iterations){
  iteration_df[,i] <- base::sample(1:length(gene_pool), 
                                   nr_recl_genes, 
                                   replace = FALSE)
  colnames(iteration_df)[i] <- i
}
set.seed(37)

print(iteration_df[1:10,1:2])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ITERATIONS

print("nr_cores")
print(nr_cores)
print("iterations")
print(iterations)
print("cut_off_prop")
print(cut_off_prop)
print("k_graph")
print(k_graph)
print("resolution_louvain")
print(resolution_louvain)

# run re-clustering and score function at once for i iterations
res_df_list <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = permuting_reclustering_scores_sce,
  sce = sce_sub, # use subsetted SCE without subclustering genes
  iteration_df = iteration_df,
  resolution_louvain = resolution_louvain,
  k_graph = k_graph,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

# for testing without mclapply
# res_df_list <- lapply(
#   X = as.list(c(1:iterations)),
#   FUN = permuting_reclustering_scores_sce,
#   sce = sce_sub,
#   iteration_df = iteration_df,
#   k_graph = k_graph,
#   resolution_louvain = resolution_louvain)

# there can be a warning about tied ks for SNNgraph, but with a set seed, this 
# should be fine 
# (https://rdrr.io/github/LTLA/BiocNeighbors/man/BiocNeighbors-ties.html)

score_df <- dplyr::bind_rows(res_df_list)
score_df$cut_off_prop <- base::rep(cut_off_prop, nrow(score_df))
print(head(score_df))

#-------------------------------------------------------------------------------

base::saveRDS(score_df, snakemake@output[["perm_score_df"]])

utils::sessionInfo()
