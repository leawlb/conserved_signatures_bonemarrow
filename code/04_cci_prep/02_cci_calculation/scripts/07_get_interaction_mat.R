library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
perc_df <- readRDS(file = snakemake@input[["perc_df"]])
lrdb <- readRDS(file = snakemake@input[["lrdb"]])
min_perc <- snakemake@params[["min_perc"]]

# use main pipeline 
source(snakemake@params[["main_functions"]])
# user helper functions
source(snakemake@params[["help_functions"]])

#-------------------------------------------------------------------------------
# helper function for quick and easy preparation of datasheet and counts matrix
prep_list <- prepare_extraction(sce = sce, assay_use = "downsampled")

idents_list <- as.list(prep_list$idents)

# helper function for implementing gene cutoff according to percentage of cells
expr_genes_list <- lapply(X = idents_list, 
                          perc_df = perc_df, 
                          cutoff = min_perc, 
                          FUN = perc_expr_genes)
    
names(expr_genes_list) <- prep_list$idents

# use prepared data to extract interaction matrix
# this function was adjusted from Adrien's function interactionmatrix()

# main function from main pipeline
interaction_mat <- extract_matrix(counts = prep_list$countsmatrix,
                                  datasheet = prep_list$datasheet,
                                  expr_genes = expr_genes_list,
                                  interactions = lrdb
                                  )

table(is.na(interaction_mat))

saveRDS(interaction_mat, snakemake@output[["interaction_mat"]])
saveRDS(prep_list$datasheet, snakemake@output[["datasheet"]])