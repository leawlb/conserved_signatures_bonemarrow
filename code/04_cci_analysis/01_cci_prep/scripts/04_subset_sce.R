library(SingleCellExperiment)

#-------------------------------------------------------------------------------

lrdb <- readRDS(file = snakemake@input[["lrdb_input"]])
sce <- readRDS(file = snakemake@input[["sce_input"]])

# subset by gene IDs
ids_lrdb <- unique(c(lrdb$ligand_ensembl_gene_id, 
                     lrdb$receptor_ensembl_gene_id))
print(length(ids_lrdb))
table(rowData(sce)$ENSMUS_ID %in% ids_lrdb)

print(nrow(sce))
sce <- sce[rowData(sce)$ENSMUS_ID %in% ids_lrdb,] 
print(nrow(sce))

saveRDS(sce, file = snakemake@output[["sce_output"]])