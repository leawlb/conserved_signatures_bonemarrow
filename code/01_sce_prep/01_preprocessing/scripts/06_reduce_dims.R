#-------------------------------------------------------------------------------
# authors: Amy Danson, Lea Wölbert
# extract highly variable genes and reduce dimensions for QC purposes

set.seed(37)
source(file = snakemake@params[["functions"]])

#-------------------------------------------------------------------------------
sce <- readRDS(file = snakemake@input[["sce_input"]])
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# get highly variably genes
# use own function summarizing the dimensionality reduction steps
sce <- reduce_dims(sce, nr_hvgs = nr_hvgs)

saveRDS(sce, file = snakemake@output[["sce_output"]])

sessionInfo()