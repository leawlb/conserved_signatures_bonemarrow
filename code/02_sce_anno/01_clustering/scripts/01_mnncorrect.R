#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(batchelor, quietly = TRUE)
source(file = snakemake@params[["functions"]])
set.seed(37)

#-------------------------------------------------------------------------------

# load object and prepare batch vector for correctExperiments() + preparation
sce <- readRDS(file = snakemake@input[["sce_input"]])

# load parameters from config
# which colData to use as batch
batch_use <- snakemake@params[["batch_use"]]
batch_pos <- which(colnames(colData(sce)) == batch_use)
# how many hvgs for batch correction
nr_hvgs_batch_correction <- snakemake@params[["nr_hvgs_batch_correction"]] 
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# rename or copy the already existant Reduceddims and assays for testing
reducedDimNames(sce)[grep("PCA", reducedDimNames(sce))] <- "PCA_before_BC"
reducedDimNames(sce)[grep("UMAP", reducedDimNames(sce))] <- "UMAP_before_BC"
assays(sce)$counts_before_BC <- assays(sce)$counts
assays(sce)$logcounts_before_BC <- assays(sce)$logcounts

# add info on batch to SCE
sce$Batch_type_used <- rep(batch_use, ncol(sce))

#-------------------------------------------------------------------------------

# fastMNN recommends use of MultiBatchNorm, output in logcounts
renorm_sce <- multiBatchNorm(sce, batch = colData(sce)[,batch_pos])
stopifnot(colnames(renorm_sce) == colnames(sce))
assays(sce)$logcounts <- assays(renorm_sce)$logcounts

print(sce)

#-------------------------------------------------------------------------------

# calculate hvgs 
gene_var <- modelGeneVar(sce)
hvgs_for_batch_correction <- getTopHVGs(gene_var, n = nr_hvgs_batch_correction)
hvgs <- getTopHVGs(gene_var, n = nr_hvgs)

#-------------------------------------------------------------------------------

print("starting MNN")

# correct using MultiBatchNorm values in logcounts assay
sce_bc <- correctExperiments(sce,
                             correct.all = TRUE,
                             PARAM = FastMnnParam(cos.norm = FALSE), # since it was previously normalized
                             batch = colData(sce)[,batch_pos],
                             subset.row = hvgs_for_batch_correction, 
                             assay.type = "logcounts") 
sce_bc$Correction_method <- rep("FastMNN", ncol(sce_bc))

print("done MNN")
print(sce_bc)

#-------------------------------------------------------------------------------
  
# rename it "PCA" because orig PCA is called PCA_before 
# so all following functions use corrected PCA automatically
reducedDimNames(sce_bc)[reducedDimNames(sce_bc) == "corrected"] <- "PCA"
colnames(reducedDim(sce_bc, type = "PCA")) <- c(1:ncol(
  reducedDim(sce_bc, type = "PCA")))
colnames(reducedDim(sce_bc, type = "PCA")) <- paste0(
  "PC", colnames(reducedDim(sce_bc, type = "PCA")))
  
# calculate UMAP based on corrected PCA
# no new hvg determination because logcounts are not corrected
seeds_umap <- snakemake@params[["seeds_umap"]]
if(sce_bc$Fraction_ID[1] == "hsc"){
  seed <- seeds_umap[["hsc"]]
}else if(sce_bc$Fraction_ID[1] == "str"){
  seed <- seeds_umap[["str"]]
}
set.seed(seed)
print(seed)
  
sce_bc <- runUMAP(sce_bc, dimred = "PCA", subset_row = hvgs)

warnings()
print(sce_bc)

#-------------------------------------------------------------------------------
saveRDS(sce_bc, file = snakemake@output[["sce_output"]])