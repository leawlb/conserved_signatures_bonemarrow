set.seed(37)

#-------------------------------------------------------------------------------

cci <- readRDS(file = snakemake@input[["cci_input"]])
idi_list <- readRDS(file = snakemake@input[["idi_list"]])

# use main pipeline 
source(snakemake@params[["metrics_functions"]])

# calculate interactions per identity

# this extracts all detected ligands or receptors per identity
ident_lrs_list <- lapply(idi_list[[2]], extract_lrs_info)

# this is more or less equivalent to the overviews of the other objects
ident_nrlrs_list <- extract_lrs_nrs(cci = cci, 
                                    lrs_list = ident_lrs_list)

saveRDS(ident_lrs_list, snakemake@output[["ident_lrs_list"]])
saveRDS(ident_nrlrs_list, snakemake@output[["ident_nrlrs_list"]])