set.seed(37)

#-------------------------------------------------------------------------------

cci <- readRDS(file = snakemake@input[["cci_input"]])
sce <- readRDS(file = snakemake@input[["sce_input"]])


# use main pipeline 
source(snakemake@params[["metrics_functions"]])

# extract nr of interactions per identity pairs and more
ipi_list <- extract_ident_pair_info(cci = cci, sce = sce)

saveRDS(ipi_list, snakemake@output[["ipi_list"]])