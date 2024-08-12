#-------------------------------------------------------------------------------

# export SCE as adata.h5ad object using zellkonverter and basilisk
# use ../../envs/zellkonverter_export.yml

library(zellkonverter, quietly = TRUE)
library(basilisk, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)

set.seed(37)

sce <- base::readRDS(snakemake@input[["sce_input"]])
#sce <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10")

#-------------------------------------------------------------------------------

# keeping only certain cols
colData(sce) <- colData(sce)[,which(colnames(colData(sce)) %in% 
                               c("batch",
                                 "Barcode", 
                                 "Object_ID", 
                                 "Species_ID",
                                 "Age_ID",
                                 "Fraction_ID",
                                 "Antibody_combination",
                                 "celltypes"))]
print(head(colData(sce)))

# keeping only raw counts and logcounts before batch correction (BC)
# these logcounts are normalised (multibatchnorm)
assays(sce) <- list(SummarizedExperiment::assays(sce)[["counts_before_BC"]], 
                    SummarizedExperiment::assays(sce)[["logcounts"]])
names(SummarizedExperiment::assays(sce)) <- c("counts", "logcounts")

# keeping only PCA (batch corrected)
SingleCellExperiment::reducedDims(sce) <- list(
  SingleCellExperiment::reducedDims(sce)[["PCA"]])
names(SingleCellExperiment::reducedDims(sce)) <- c("PCA")

# completely remove metadata to circumvent a bug
S4Vectors::metadata(sce) <- list()

# remove rotations from rowData
rowData(sce) <- rowData(sce)[, c(2, 3, 4)]
print(head(rowData(sce)))

#-------------------------------------------------------------------------------

# save and convert in one step
zellkonverter::writeH5AD(sce,
                         snakemake@output[["adata_output"]], 
                         X_name = "counts", 
                         verbose = TRUE)

utils::sessionInfo()
