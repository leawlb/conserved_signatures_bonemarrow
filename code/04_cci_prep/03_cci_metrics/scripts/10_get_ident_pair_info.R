set.seed(37)

#-------------------------------------------------------------------------------

interaction_list <- readRDS(file = snakemake@input[["interaction_list"]])

age_curr <- snakemake@wildcards[["age"]]
print(age_curr)
species_curr <- snakemake@wildcards[["species"]]
print(species_curr)

# use main pipeline 
source(snakemake@params[["main_functions"]])

# extract nr of interactions per identity pairs and more
ident_pair_info <- extract_ident_pair_info(int_list = interaction_list)

ident_pair_info$overview$age <- age_curr
ident_pair_info$overview$species <- species_curr
ident_pair_info$overview$condition <- "main"

saveRDS(ident_pair_info, snakemake@output[["ident_pair_info"]])