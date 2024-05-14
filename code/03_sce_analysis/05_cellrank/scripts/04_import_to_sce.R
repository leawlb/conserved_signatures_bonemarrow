# convert HSC adata back to SCE object

#-------------------------------------------------------------------------------
library(zellkonverter)
library(basilisk)

#-------------------------------------------------------------------------------

# convert and load in one step from zellkonverter
sce <- readH5AD(snakemake@input[["adata_input"]])

# save SCE
print(snakemake@output[["sce_output"]])
saveRDS(sce, snakemake@output[["sce_output"]])