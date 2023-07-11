
#-------------------------------------------------------------------------------

ident_pair_info <- readRDS(file = snakemake@input[["ident_pair_info"]])

age_curr <- snakemake@wildcards[["age"]]
print(age_curr)
species_curr <- snakemake@wildcards[["species"]]
print(species_curr)

# use main pipeline 
source("../../source/cci_functions_calculation_main.R")

# calculate interactions per identity
ident_info <- extract_ident_info(ipi_list = ident_pair_info)

ident_info$overview$age <- age_curr
ident_info$overview$species <- species_curr
ident_info$overview$condition <- "main"

saveRDS(ident_info, snakemake@output[["ident_info"]])