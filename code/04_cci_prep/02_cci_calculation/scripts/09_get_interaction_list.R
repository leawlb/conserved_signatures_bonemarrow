set.seed(37)

#-------------------------------------------------------------------------------

interaction_score_list <- readRDS(file = snakemake@input[["interaction_scoring"]])
lrdb <- readRDS(file = snakemake@input[["lrdb"]])

age_curr <- snakemake@wildcards[["age"]]
species_curr <- snakemake@wildcards[["species"]]

# use main pipeline 
source(snakemake@params[["main_functions"]])

interaction_list <- make_interaction_list(
  interaction_scores = interaction_score_list,
  lrdb = lrdb)

interaction_list$Identities$Age <- rep(age_curr, 
                                       nrow(interaction_list$Identities))
interaction_list$Identities$Species <- rep(species_curr, 
                                           nrow(interaction_list$Identities))

saveRDS(interaction_list, snakemake@output[["interaction_list"]])