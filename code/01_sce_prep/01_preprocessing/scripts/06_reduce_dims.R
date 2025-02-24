#-------------------------------------------------------------------------------
# authors: Amy Danson, Lea Wölbert
# extract highly variable genes and reduce dimensions for QC purposes

RNGkind("L'Ecuyer-CMRG")
set.seed(37)


source(file = snakemake@params[["functions"]])

#-------------------------------------------------------------------------------
sce <- base::readRDS(file = snakemake@input[["sce_input"]])
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# use own function summarizing the basic dimensionality reduction steps
sce <- reduce_dims(sce, nr_hvgs = nr_hvgs)

base::saveRDS(sce, file = snakemake@output[["sce_output"]])

utils::sessionInfo()