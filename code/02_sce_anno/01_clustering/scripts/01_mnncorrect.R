#-------------------------------------------------------------------------------
# MNN batch correction (fastMNN)

library(scran, quietly = TRUE)
library(scater, quietly = TRUE)
library(batchelor, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])

batch_use <- snakemake@params[["batch_use"]]
batch_pos <- which(colnames(colData(sce)) == batch_use)

nr_hvgs_BC <- snakemake@params[["nr_hvgs_BC"]] 
nr_hvgs <- snakemake@params[["nr_hvgs"]]

#-------------------------------------------------------------------------------
# rename the already existent reduced dims and assays (to keep them)
reducedDimNames(sce)[grep("PCA", reducedDimNames(sce))] <- "PCA_before_BC"
reducedDimNames(sce)[grep("UMAP", reducedDimNames(sce))] <- "UMAP_before_BC"
assays(sce)$counts_before_BC <- assays(sce)$counts
assays(sce)$logcounts_before_BC <- assays(sce)$logcounts

# add info on batch type used
sce$Batch_type_used <- rep(batch_use, ncol(sce))

#-------------------------------------------------------------------------------

# fastMNN recommends use of MultiBatchNorm
# output in logcounts and logcounts_batchnorm (for reference)
renorm_sce <- batchelor::multiBatchNorm(sce, batch = colData(sce)[,batch_pos])
stopifnot(colnames(renorm_sce) == colnames(sce))
assays(sce)$logcounts_batchnorm <- assays(renorm_sce)$logcounts
assays(sce)$logcounts <- assays(renorm_sce)$logcounts

print(sce)

#-------------------------------------------------------------------------------

# calculate hvgs 
gene_var <- scran::modelGeneVar(sce)
hvgs_BC <- scran::getTopHVGs(gene_var, n = nr_hvgs_BC)
hvgs <- scran::getTopHVGs(gene_var, n = nr_hvgs)

#-------------------------------------------------------------------------------

print("starting MNN")

# correct using MultiBatchNorm values in logcounts assay
sce_bc <- batchelor::correctExperiments(
  sce,
  correct.all = TRUE,
  PARAM = FastMnnParam(cos.norm = FALSE), # since it was previously normalized
  batch = colData(sce)[,batch_pos],
  subset.row = hvgs_BC, 
  assay.type = "logcounts") 
sce_bc$Correction_method <- rep("FastMNN", ncol(sce_bc))

print("done MNN")
print(sce_bc)
warnings() # warnings about "useNames = NA is deprecated" which seems ignorable

#-------------------------------------------------------------------------------
  
# rename it "PCA" because orig PCA is now called PCA_before 
# so all following functions use corrected PCA automatically
reducedDimNames(sce_bc)[reducedDimNames(sce_bc) == "corrected"] <- "PCA"
colnames(reducedDim(sce_bc, type = "PCA")) <- c(1:ncol(
  reducedDim(sce_bc, type = "PCA")))
colnames(reducedDim(sce_bc, type = "PCA")) <- paste0(
  "PC", colnames(reducedDim(sce_bc, type = "PCA")))
  
# calculate UMAP based on corrected PCA
# no new hvg determination because logcounts_batchnorm is not corrected anyway
seeds_umap <- snakemake@params[["seeds_umap"]]
if(sce_bc$Fraction_ID[1] == "hsc"){
  seed <- seeds_umap[["hsc"]]
}else if(sce_bc$Fraction_ID[1] == "str"){
  seed <- seeds_umap[["str"]]
}

set.seed(seed)
print(seed)
sce_bc <- scater::runUMAP(sce_bc, dimred = "PCA", subset_row = hvgs)
print(sce_bc)

#-------------------------------------------------------------------------------
saveRDS(sce_bc, file = snakemake@output[["sce_output"]])

sessionInfo()