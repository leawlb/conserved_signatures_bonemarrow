#-------------------------------------------------------------------------------
# MNN batch correction (fastMNN)

library(scran, quietly = TRUE)
library(scater, quietly = TRUE)
library(batchelor, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# objects and params

sce <- base::readRDS(file = snakemake@input[["sce_input"]])

batch_use <- snakemake@params[["batch_use"]]
batch_pos <- which(colnames(colData(sce)) == batch_use)

nr_hvgs_BC <- snakemake@params[["nr_hvgs_BC"]] 
nr_hvgs <- snakemake@params[["nr_hvgs"]]

seeds_umap <- snakemake@params[["seeds_umap"]]
fraction_curr <- snakemake@wildcards[["fraction"]]

print(seeds_umap)
print(fraction_curr)
print(nr_hvgs_BC)

seed <- seeds_umap[[fraction_curr]]

print(seed)

#-------------------------------------------------------------------------------
# rename the already existent reduced dims and assays (to keep them)
SingleCellExperiment::reducedDimNames(sce)[base::grep(
  "PCA", SingleCellExperiment::reducedDimNames(sce))] <- "PCA_before_BC"
SingleCellExperiment::reducedDimNames(sce)[base::grep(
  "UMAP", SingleCellExperiment::reducedDimNames(sce))] <- "UMAP_before_BC"

SummarizedExperiment::assays(sce)$counts_before_BC <- SummarizedExperiment::assays(sce)$counts
SummarizedExperiment::assays(sce)$logcounts_before_BC <- SummarizedExperiment::assays(sce)$logcounts

# add info on batch type used
sce$Batch_type_used <- base::rep(batch_use, ncol(sce))

#-------------------------------------------------------------------------------

# fastMNN recommends use of MultiBatchNorm
# output in "logcounts" for downstream use and
# "logcounts_batchnorm" for reference
set.seed(37)
renorm_sce <- batchelor::multiBatchNorm(sce, 
                                        batch = colData(sce)[,batch_pos])

SummarizedExperiment::assays(sce)$logcounts_batchnorm <- SummarizedExperiment::assays(renorm_sce)$logcounts
SummarizedExperiment::assays(sce)$logcounts <- SummarizedExperiment::assays(renorm_sce)$logcounts

print(SingleCellExperiment::reducedDimNames(sce))
print(SummarizedExperiment::assays(sce))

#-------------------------------------------------------------------------------

# calculate hvgs for batch correction and UMAP
gene_var <- scran::modelGeneVar(sce)
hvgs_BC <- scran::getTopHVGs(gene_var, n = nr_hvgs_BC)
hvgs <- scran::getTopHVGs(gene_var, n = nr_hvgs)

#-------------------------------------------------------------------------------

print("starting MNN")
set.seed(37)

# correct using MultiBatchNorm values in logcounts assay
sce_bc <- batchelor::correctExperiments(
  sce,
  correct.all = TRUE,
  PARAM = FastMnnParam(
    cos.norm = FALSE, # since it is already normalized
    k=20, # default
    d=50), # default
  batch = colData(sce)[,batch_pos],
  subset.row = hvgs_BC,  
  assay.type = "logcounts") 

print("done MNN")
print(sce_bc)
warnings() # warnings about "useNames = NA is deprecated" which seems ignorable

#-------------------------------------------------------------------------------

# store corrected values in "PCA" slot which will be used by downstream
# original PCA is still stored in "PCA_before_BC" slot
SingleCellExperiment::reducedDimNames(sce_bc)[
  SingleCellExperiment::reducedDimNames(sce_bc) == "corrected"] <- "PCA"

# rename PCA colnames
colnames(SingleCellExperiment::reducedDim(sce_bc, type = "PCA")) <- c(1:ncol(
  SingleCellExperiment::reducedDim(sce_bc, type = "PCA")))
colnames(SingleCellExperiment::reducedDim(sce_bc, type = "PCA")) <- base::paste0(
  "PC", colnames(SingleCellExperiment::reducedDim(sce_bc, type = "PCA")))

print(SingleCellExperiment::reducedDim(sce_bc, type = "PCA")[1:5, 1:5])

#-------------------------------------------------------------------------------
# run umap with fraction-specific seed for nice looking plots
# hvgs are calculated from normalised logcounts
set.seed(seed)
sce_bc <- scater::runUMAP(sce_bc,
                          dimred = "PCA", 
                          subset_row = hvgs)

set.seed(37)

print(head(reducedDim(sce_bc, type = "UMAP")))

#-------------------------------------------------------------------------------

base::saveRDS(sce_bc, file = snakemake@output[["sce_output"]])

utils::sessionInfo()