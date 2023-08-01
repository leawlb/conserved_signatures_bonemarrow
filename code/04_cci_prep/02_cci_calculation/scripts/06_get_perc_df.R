library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])

# use main pipeline 
source(snakemake@params[["main_functions"]])

# use downsampled for main pipeline
exp_gene_perc_df <- expr_cells_perc(sce, assay_use = "downsampled")

print(head(exp_gene_perc_df))

saveRDS(exp_gene_perc_df, snakemake@output[["perc_df"]])