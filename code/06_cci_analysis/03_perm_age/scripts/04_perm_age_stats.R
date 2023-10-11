library(tidyverse)
library(parallel)
library(data.table)

set.seed(37)

source(snakemake@params[["perm_functions"]])

cci_yng <- readRDS(file = snakemake@input[["cci_yng"]])
cci_old <- readRDS(file = snakemake@input[["cci_old"]])

perm_scores_yng <- readRDS(file = snakemake@input[["perm_scores_yng"]])
perm_scores_old <- readRDS(file = snakemake@input[["perm_scores_old"]])
#perm_scores_yng <- perm_score_df_lists$yng # for manual run
#perm_scores_old <- perm_score_df_lists$old

nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]
min_sp_age_perc <- snakemake@params[["min_sp_age_perc"]]

# min_sp_age_perc <- 0.1

#-------------------------------------------------------------------------------
# sort lists so that the items contain their own names
# important for interaction_stats_age()

perm_scores_yng_named <- list()
for(i in names(perm_scores_yng)){
  perm_scores_yng_named[[i]]$name <- i
  perm_scores_yng_named[[i]]$df <- perm_scores_yng[[i]]
}

perm_scores_old_named <- list()
for(i in names(perm_scores_old)){
  perm_scores_old_named[[i]]$name <- i
  perm_scores_old_named[[i]]$df <- perm_scores_old[[i]]
}

#-------------------------------------------------------------------------------
# iterate through cell type pairs (ctps) 

ctps <- as.list(names(perm_scores_yng))
# only possible for cell type pairs that were detected in both
ctps <- ctps[ctps %in% names(perm_scores_old)]
ctps <- ctps[ctps %in% colnames(cci_yng$Score)]
ctps <- ctps[ctps %in% colnames(cci_old$Score)]

# get statistics for permutation tests for each interaction in each ctp 
# statistics of score deltas = score old - score yng
res_lists <- mclapply(
  X = ctps,
  min_sp_age_perc = min_sp_age_perc,
  perm_score_list_yng = perm_scores_yng_named, 
  perm_score_list_old = perm_scores_old_named,
  cci_yng = cci_yng,
  cci_old = cci_old,
  iterations = iterations,
  FUN = stats_perm_age,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

#res_lists <- lapply(
#  X = ctps, 
#  min_sp_age_perc = min_sp_age_perc,
#  perm_score_list_yng = perm_scores_yng_named, 
#  perm_score_list_old = perm_scores_old_named,
#  cci_yng = cci_yng, 
#  cci_old = cci_old,
#  iterations = iterations,
#  FUN = stats_perm_age)

warnings()

names(res_lists) <- names(perm_scores_yng)
print(names(res_lists))
print(names(res_lists[[1]]))

saveRDS(res_lists, snakemake@output[["res_lists"]])