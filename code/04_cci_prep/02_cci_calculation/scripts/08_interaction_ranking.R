library(dplyr)
set.seed(37)

#-------------------------------------------------------------------------------

interaction_mat <- readRDS(file = snakemake@input[["interaction_mat"]])
datasheet <- readRDS(file = snakemake@input[["datasheet"]])
top_level <- snakemake@params[["top_level"]]

# use main pipeline 
source(snakemake@params[["main_functions"]])

# this function was adjusted from Adrien's interactionranking() function
interaction_rank_list <- interaction_ranking(interaction_mat = interaction_mat, 
                                             datasheet = datasheet,
                                             top_level = top_level)

saveRDS(interaction_rank_list, snakemake@output[["interaction_ranking"]])