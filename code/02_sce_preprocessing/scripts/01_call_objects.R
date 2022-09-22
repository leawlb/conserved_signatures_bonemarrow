#-------------------------------------------------------------------------------
# author: Amy Danson, Lea Wölbert
# date: 2022-09-19
# call objects from raw cellranger output matrices and add metadata

# load objects and libraries
library(SingleCellExperiment, quietly = TRUE) 
metadata <- read.csv(file = "../metadata.csv", header = TRUE, sep = ";")

# call SCE objects from the cellranger output matrices
sce <- read10xCounts(file = snakemake@input[["cr_output"]], col.names = TRUE, 
                     type = "sparse" )

# reconstruct the object ID of the currently loaded object from the wildcards 
# they should be identical by design
object_id_curr <- paste(
  snakemake@wildcards[["Species"]],
  snakemake@wildcards[["Age"]],
  snakemake@wildcards[["Fraction"]],
  snakemake@wildcards[["Sample"]],
  sep = "_"
)

# retrieve only information regarding the loaded object
metadata_curr <- metadata[metadata$object_id == object_id_curr,]

# add cell-level metadata to each barcode 
for(i in colnames(metadata_curr)){
  colData(sce)[i] <- rep(metadata_curr[,i], ncol(sce))
}

# change rownames to Symbols for better readability
rownames(sce) <- rowData(sce)$Symbol

# save as new object
saveRDS(sce, file = snakemake@output[["sce_01"]])

# testing purposes
#sce_01 <- read10xCounts("/omics/odcf/analysis/OE0538_projects/DO-0008/mmus/mmus_old/croutput_files/mmus_OLD1_bm2/outs/raw_feature_bc_matrix", col.names = TRUE, type = "sparse" )
#object_id_curr <- "mmus_old_hsc_1.0"
name_curr <- "mmus_old_hsc_1.0"
#sce <- sce_01
#saveRDS(sce, file = "../data/preprocessing/call/sce_mmus_old_hsc_1.0_01")
sce_02 <- readRDS(file = "../data/preprocessing/drop/sce_mmus_old_hsc_1.0_02")
#setwd("~/Interspecies_BM/preprocessing")  
  
  