#-------------------------------------------------------------------------------

# export SCE as adata.h5ad object using zellkonverter and basilisk

library(zellkonverter)
library(basilisk)
library(SingleCellExperiment)

sce <- readRDS(snakemake@input[["sce_input"]])
#sce <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10")

#-------------------------------------------------------------------------------

# keeping only certain cols
colData(sce) <- colData(sce)[,c(1, 3, 5, 8, 10, 13, 14, 37, 42)]

# keeping only raw counts and logcounts before batch correction (BC)
assays(sce) <- list(assays(sce)[["counts_before_BC"]], 
                    assays(sce)[["logcounts_before_BC"]])
names(assays(sce)) <- c("counts", "logcounts")

# keeping only PCA (batch corrected)
reducedDims(sce) <- list(reducedDims(sce)[["PCA"]])
names(reducedDims(sce)) <- c("PCA")

# completely remove metadata to circumvent a bug
metadata(sce) <- list()

# remove rotations from rowData
rowData(sce) <- rowData(sce)[, c(2, 3, 4)]

#-------------------------------------------------------------------------------

# save and convert in one step
zellkonverter::writeH5AD(sce, snakemake@output[["adata_output"]], 
                         X_name = "counts", verbose = TRUE)
