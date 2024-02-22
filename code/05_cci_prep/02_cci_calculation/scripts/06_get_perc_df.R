library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])

source(snakemake@params[["main_functions"]])

# use downsampled counts 
exp_gene_perc_df <- perc_expr_cells(sce, assay_use = "downsampled")

print(head(exp_gene_perc_df))

saveRDS(exp_gene_perc_df, snakemake@output[["perc_df"]])