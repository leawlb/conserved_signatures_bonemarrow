
#-------------------------------------------------------------------------------

interaction_rank_list <- readRDS(file = snakemake@input[["interaction_ranking"]])
lrdb <- readRDS(file = snakemake@input[["lrdb"]])

# use main pipeline 
source(snakemake@params[["main_functions"]])

interaction_list <- make_interaction_list(
  interaction_ranking = interaction_rank_list,
  lrdb = lrdb)

saveRDS(interaction_list, snakemake@output[["interaction_list"]])