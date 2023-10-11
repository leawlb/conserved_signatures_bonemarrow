library(tidyverse)
library(SingleCellExperiment)
library(parallel)
library(data.table)

RNGkind("L'Ecuyer-CMRG") # for mclapply, generation of same random numbers

set.seed(37)

source(snakemake@params[["perm_functions"]])
source(snakemake@params[["main_functions"]])

sce_yng <- readRDS(file = snakemake@input[["sce_yng_input"]])
sce_old <- readRDS(file = snakemake@input[["sce_old_input"]])
lrdb <- readRDS(file = snakemake@input[["lrdb"]])

iterations <-  snakemake@params[["iterations"]]
top_level <- snakemake@params[["top_level"]]
min_perc <- snakemake@params[["min_perc"]]
nr_cores <- snakemake@params[["nr_cores"]]

species_curr <- snakemake@wildcards[["species"]]
ages <- snakemake@params[["ages"]]

#-------------------------------------------------------------------------------
# prepare SCE objects (make them smaller)

# only keep cell types that are existent in both ages after SCE prep  
sce_yng <- sce_yng[,sce_yng$Identity %in% unique(sce_old$Identity)]
sce_old <- sce_old[,sce_old$Identity %in% unique(sce_yng$Identity)]

# remove unnecessary assays and data
assays(sce_yng) <- list("downsampled" = assays(sce_yng)$downsampled)
assays(sce_old) <- list("downsampled" = assays(sce_old)$downsampled)
reducedDims(sce_yng) <- NULL
metadata(sce_yng) <- list()
reducedDims(sce_old) <- NULL
metadata(sce_old) <- list()

colData(sce_yng) <- colData(sce_yng)[,which(colnames(colData(sce_yng)) %in% 
                                              c("Age_ID", "Species_ID",
                                                "Identity", "Assignment", 
                                                "name", "temp_ident"))]
colData(sce_old) <- colData(sce_old)[,which(colnames(colData(sce_old)) %in% 
                                              c("Age_ID", "Species_ID",
                                                "Identity", "Assignment",
                                                "name", "temp_ident"))]

# more complicated to keep the sequence of the factors the same
sce_yng$Identity <- factor(
  sce_yng$Identity, 
  levels = levels(sce_yng$Identity)[
    which(levels(sce_yng$Identity) %in% unique(sce_yng$Identity))])

sce_old$Identity <- factor(
  sce_old$Identity, 
  levels = levels(sce_old$Identity)[
    which(levels(sce_old$Identity) %in% unique(sce_old$Identity))])

# prepare lrdb (make it smaller)
lrdb <- lrdb[,colnames(lrdb) %in% c("lr_pair", "ligand_gene_symbol", 
                                    "receptor_gene_symbol", 
                                    "ligand_ensembl_gene_id",
                                    "receptor_ensembl_gene_id")]

---------------------------------------------------------------------
# permute age labels and calculate CCI scores in one step
# use mclapply for parallel processing 

#iterations <- 50
# this function returns a list of length = iterations with permuted "yng"
# and "old" CCI objects
cci_list_ages <- mclapply(
  X = c(1:iterations),
  FUN = permute_age,
  sce_old = sce_old, 
  sce_yng = sce_yng,
  lrdb = lrdb, 
  min_perc = min_perc,
  top_level = top_level, 
  mc.preschedule = TRUE, 
  mc.cores = nr_cores,
  mc.silent = TRUE)

#-------------------------------------------------------------------------------
# sort list by ages for separate processing

names(cci_list_ages) <- paste0("iteration_", c(1:iterations))

cci_list_ages_sorted <- list()
for(a in ages){
  for(i in names(cci_list_ages)){
    cci_list_ages_sorted[[a]][[i]] <- cci_list_ages[[i]][[a]]
  }
}

#-------------------------------------------------------------------------------
# convert list of n permuted CCIs into list of perm score dfs with 
# rows = interactions, cols = iterations, per interacting cell type pair (ctp)

perm_score_df_lists <- lapply(
  X = cci_list_ages_sorted, 
  FUN = convert_perm_df, 
  lrdb = lrdb)

names(perm_score_df_lists) <- names(cci_list_ages_sorted)

saveRDS(perm_score_df_lists[["yng"]], snakemake@output[["perm_scores_yng"]])
saveRDS(perm_score_df_lists[["old"]], snakemake@output[["perm_scores_old"]])

time <- Sys.time()
print("after saving")
print(time - time1)
sort(sapply(ls(),function(x){object.size(get(x))})) 
