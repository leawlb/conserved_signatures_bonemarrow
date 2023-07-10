library(SingleCellExperiment)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])

# use main pipeline 
source("../../source/cci_functions_calculation_main.R")

# use downsampled for main pipeline
exp_gene_perc_df <- expr_cells_perc(sce, assay_use = "downsampled")

print(head(exp_gene_perc_df))

saveRDS(exp_gene_perc_df, snakemake@output[["perc_df"]])