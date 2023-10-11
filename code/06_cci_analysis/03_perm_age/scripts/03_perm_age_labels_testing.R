time1 <- Sys.time()
print(time1)

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

time <- Sys.time()
print("after loading")
print(time - time1)
sort(sapply(ls(),function(x){object.size(get(x))})) 

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

# more complicate to keep the sequence of the factors the same
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


time <- Sys.time()
print("after prep") 
print(time - time1)
sort(sapply(ls(),function(x){object.size(get(x))})) 

#-------------------------------------------------------------------------------
# there seems to be an effect where a larger score delta also results in 
# larger permuted score deltas so that the entire permuted distribution can be
# shifted significantly away from 0
# this effect may be caused by large differences in cell number per identity
# between young and old SCE datasets (or something else, no idea)


# subset SCE so that each ctp will have the same number of cells (45, which is the minimum)
# if different nrs of cells is at fault, the effect should be abolished
# THIS IS FOR TESTING
#ctps <- unique(unfactor(sce_old$Identity))
#table(sce_old$Identity)
#table(sce_yng$Identity)

#for(c in ctps){
#  old_pos <- which(sce_old$Identity == c)
#  old_pos_other <- which(sce_old$Identity != c)
#  yng_pos <- which(sce_yng$Identity == c)
#  yng_pos_other <- which(sce_yng$Identity != c)
#
#  sce_old <- sce_old[, c(old_pos[1:45], old_pos_other)]
#  sce_yng <- sce_yng[, c(yng_pos[1:45], yng_pos_other)]
#}

#-------------------------------------------------------------------------------
# permute cell type labels to see if the effect is cell type specific 
# this will keep the numbers of cells the same
# but remove cell type-specific effect
# if different nrs of cells are at fault, the effect should persist
# THIS IS FOR TESTING
#perm_cts_old <- sample(c(1:ncol(sce_old)), size = ncol(sce_old))
#perm_cts_yng <- sample(c(1:ncol(sce_yng)), size = ncol(sce_yng))

#sce_old$Identity <- sce_old$Identity[perm_cts_old]
#sce_yng$Identity <- sce_yng$Identity[perm_cts_yng]

#-------------------------------------------------------------------------------
## check if it works when I remove the cell types below 100
# not really

#print(table(sce_old$Identity))
#print(table(sce_yng$Identity))

#ct_rm <- names(table(sce_old$Identity))[which(table(sce_old$Identity) < 100)]
#ct_rm <- c(ct_rm, names(table(sce_yng$Identity))[which(table(sce_yng$Identity) < 100)])
#ct_rm

#sce_yng <- sce_yng[,-which(sce_yng$Identity %in% ct_rm)]
#sce_old <- sce_old[,-which(sce_old$Identity %in% ct_rm)]

## TODO: change downsampling to downsampling of all cell types to the lowest
## cell type and see what happens :(


#-------------------------------------------------------------------------------

sce_yng$Identity <- factor(
  sce_yng$Identity, 
  levels = levels(sce_yng$Identity)[
    which(levels(sce_yng$Identity) %in% unique(sce_yng$Identity))])

sce_old$Identity <- factor(
  sce_old$Identity, 
  levels = levels(sce_old$Identity)[
    which(levels(sce_old$Identity) %in% unique(sce_old$Identity))])

print(table(sce_old$Identity))
print(table(sce_yng$Identity))
#-------------------------------------------------------------------------------

# check if the trend is abolished with stricter cutoff 
#min_perc <- 0.05

#-------------------------------------------------------------------------------

#print(sce_yng$Identity[1:10])
#print(sce_old$Identity[1:10])

#-------------------------------------------------------------------------------
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


#cci_list_ages <- lapply(
#  X = c(1:iterations),
#    FUN = permute_age,
#  sce_old = sce_old, 
#  sce_yng = sce_yng,
#  lrdb = lrdb, 
#  min_perc = min_perc,
#  top_level = min_perc)

time <- Sys.time()
print("after cci calc")
print(time - time1)
sort(sapply(ls(),function(x){object.size(get(x))})) 


#-------------------------------------------------------------------------------
# sort list by ages for separate processing

names(cci_list_ages) <- paste0("iteration_", c(1:iterations))

cci_list_ages_sorted <- list()
for(a in ages){
  for(i in names(cci_list_ages)){
    cci_list_ages_sorted[[a]][[i]] <- cci_list_ages[[i]][[a]]
  }
}
remove(cci_list_ages)

print(length(cci_list_ages_sorted[[1]]))
time <- Sys.time()
print("after sorting")
print(time - time1)
sort(sapply(ls(),function(x){object.size(get(x))})) 

#-------------------------------------------------------------------------------
# convert list of n permuted CCIs into list of perm score dfs with 
# rows = interactions, cols = iterations, per interacting cell type pair (ctp)

perm_score_df_lists <- lapply(
  X = cci_list_ages_sorted, 
  FUN = convert_perm_df, 
  lrdb = lrdb)

time <- Sys.time()
print("after converting")
print(time - time1)

names(perm_score_df_lists) <- names(cci_list_ages_sorted)

saveRDS(perm_score_df_lists[["yng"]], snakemake@output[["perm_scores_yng"]])
saveRDS(perm_score_df_lists[["old"]], snakemake@output[["perm_scores_old"]])

time <- Sys.time()
print("after saving")
print(time - time1)
sort(sapply(ls(),function(x){object.size(get(x))})) 
