set.seed(37)

#-------------------------------------------------------------------------------

interaction_list <- readRDS(file = snakemake@input[["interaction_list"]])
ident_pair_info <- readRDS(file = snakemake@input[["ident_pair_info"]])

age_curr <- snakemake@wildcards[["age"]]
print(age_curr)
species_curr <- snakemake@wildcards[["species"]]
print(species_curr)

# use main pipeline 
source(snakemake@params[["main_functions"]])

# calculate interactions per identity

# this extracts all detected ligands or receptors per cell type
ident_lrs_info <- lapply(ident_pair_info[[2]], extract_lrs_info)

# this is more or less equivalent to the overviews of the other objects
ident_nrlrs_info <- extract_lrs_nrs(int_list = interaction_list, 
                                    lrs_list = ident_lrs_info)
ident_nrlrs_info

ident_nrlrs_info$age <- age_curr
ident_nrlrs_info$species <- species_curr
ident_nrlrs_info$condition <- "main"

saveRDS(ident_lrs_info, snakemake@output[["ident_lrs_info"]])
saveRDS(ident_nrlrs_info, snakemake@output[["ident_nrlrs_info"]])