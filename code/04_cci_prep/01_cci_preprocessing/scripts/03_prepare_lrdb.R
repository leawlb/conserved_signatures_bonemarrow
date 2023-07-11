library(SingleCellExperiment)

#-------------------------------------------------------------------------------

lrdb <- readRDS(file = snakemake@input[["lrdb_input"]])
sce_input_path <- snakemake@input[["sce_input_path"]]

sce_list <- list()

for(i in 1:length(sce_input_path)){
  sce_list[[i]] <- readRDS(file = sce_input_path[[i]])
}

# subset by gene IDs
ids_sce <- vector()
for(i in 1:length(sce_list)){
  ids_sce <- unique(c(ids_sce, rowData(sce_list[[i]])$ENSMUS_ID))
}

print(length(ids_sce))

table(lrdb$ligand_ensembl_gene_id %in% ids_sce)
table(lrdb$receptor_ensembl_gene_id %in% ids_sce)

print(nrow(lrdb))
lrdb <- lrdb[lrdb$ligand_ensembl_gene_id %in% ids_sce,] 
lrdb <- lrdb[lrdb$receptor_ensembl_gene_id %in% ids_sce,] 
print(nrow(lrdb))

saveRDS(lrdb, file = snakemake@output[["lrdb_output"]])