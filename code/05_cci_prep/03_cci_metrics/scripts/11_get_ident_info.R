set.seed(37)

#-------------------------------------------------------------------------------

ipi_list <- readRDS(file = snakemake@input[["ipi_list"]])

# use main pipeline 
source(snakemake@params[["metrics_functions"]])

# calculate interactions per identity
idi_list <- extract_ident_info(ipi_list = ipi_list)

saveRDS(idi_list, snakemake@output[["idi_list"]])