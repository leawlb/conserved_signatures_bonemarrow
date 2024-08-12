# convert HSC adata back to SCE object

# use ../../envs/zellkonverter_import.yml

#-------------------------------------------------------------------------------
library(zellkonverter, quietly = TRUE)
library(basilisk, quietly = TRUE)

set.seed(37)

#-------------------------------------------------------------------------------

# convert and load in one step from zellkonverter
sce <- zellkonverter::readH5AD(snakemake@input[["adata_input"]])

# save SCE
print(snakemake@output[["sce_output"]])
base::saveRDS(sce, snakemake@output[["sce_output"]])