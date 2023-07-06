#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(batchelor, quietly = TRUE)
source(file = "../../source/sce_functions.R")
set.seed(37)

#-------------------------------------------------------------------------------

# load object and prepare batch vector for correctExperiments() + preparation
sce <- readRDS(file = snakemake@input[["sce_input"]])

# load parameters from config
batch_use <- snakemake@params[["batch_use"]] # which colData to use as batch
batch_pos <- which(colnames(colData(sce)) == batch_use)

nr_hvgs <- snakemake@params[["nr_hvgs"]] # how many hvgs for UMAP calc or PCA&UMAP calc
hvgs_for_batch_correction <- snakemake@params[["hvgs_for_batch_correction"]] # use hvgs?
nr_hvgs_batch_correction <- snakemake@params[["nr_hvgs_batch_correction"]] # how many hvgs for batch correction

# rename or copy the already existent Reduceddims and assays for testing
reducedDimNames(sce)[grep("PCA", reducedDimNames(sce))] <- "PCA_before_BC"
reducedDimNames(sce)[grep("UMAP", reducedDimNames(sce))] <- "UMAP_before_BC"
assays(sce)$counts_before_BC <- assays(sce)$counts
assays(sce)$logcounts_before_BC <- assays(sce)$logcounts

# determine which assay to use 
assay_use <- "logcounts"

# determine witch method to use
mnn_fast <- snakemake@params[["mnn_fast"]]

# add info on batch to SCE
sce$Batch_type_used <- rep(batch_use, ncol(sce))

#-------------------------------------------------------------------------------

if(mnn_fast == TRUE){
  
  # fastMNN recommends use of MultiBatchNorm
  # rescales size factors, output in logcounts
  renorm_sce <- multiBatchNorm(sce, batch = colData(sce)[,batch_pos])
  stopifnot(colnames(renorm_sce) == colnames(sce))
  assays(sce)$logcounts <- assays(renorm_sce)$logcounts
}

print(sce)

#-------------------------------------------------------------------------------

# calculate hvgs if required, otherwise keep NULL
if(hvgs_for_batch_correction == TRUE){
  hvgs_for_batch_correction <- modelGeneVar(sce)
  hvgs_for_batch_correction <- getTopHVGs(hvgs_for_batch_correction, 
                                          n=nr_hvgs_batch_correction)
}else{
  hvgs_for_batch_correction <- NULL  
  print("no hvgs")
}

#-------------------------------------------------------------------------------

print("starting MNN")

# run MNNfunction depending on specified method
if(mnn_fast == TRUE){
  
  print(sce)
  # correct using MultiBatchNorm values in logcounts assay
  sce_bc <- correctExperiments(sce,
                               correct.all = TRUE,
                               PARAM = FastMnnParam(cos.norm = FALSE), # since it was previously normalized
                               batch = colData(sce)[,batch_pos],
                               subset.row = hvgs_for_batch_correction, 
                               assay.type = assay_use) 
  sce_bc$Correction_method <- rep("FastMNN", ncol(sce_bc))
  # output = "corrected" (PCA) + all already existing slots
  print("done MNN")
  print(sce_bc)
  
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
  
  sce_bc <- runUMAP(sce_bc, dimred = "PCA", 
                    subset_row = hvgs_for_batch_correction)
  
}else if(mnn_fast == FALSE){
  
  # does not require MultiBatchNorm, uses its own normalization method
  sce_bc <- correctExperiments(sce, 
                               correct.all = TRUE,
                               batch = colData(sce)[,batch_pos],
                               PARAM = ClassicMnnParam(), 
                               subset.row = hvgs_for_batch_correction,
                               assay.type = assay_use)
  sce_bc$Correction_method <- rep("Classic_MNN", ncol(sce_bc))
  # output in first assay (counts)
  
  # Reduce dimensions
  seeds_umap <- snakemake@params[["seeds_umap"]]
  if(sce_bc$Fraction_ID == "hsc"){
    seed <- seeds_umap[["hsc"]]
  }else if(sce_bc$Fraction_ID == "str"){
    seed <- seeds_umap[["str"]]
  }
  set.seed(seed)
  print(seed)
  sce_bc <- reduce_dims(sce_bc, nr_hvgs = nr_hvgs)
}

warnings()

print(sce_bc)
saveRDS(sce_bc, file = snakemake@output[["sce_output"]])