#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(batchelor, quietly = TRUE)
source(file = snakemake@params[["sce_functions"]])
set.seed(37)

#-------------------------------------------------------------------------------

# load object and prepare batch vector for correctExperiments() + preparation
sce <- readRDS(file = snakemake@input[[1]])

batch_use <- snakemake@params[["batch_use"]] # which batch?
batch_pos <- which(colnames(colData(sce)) == batch_use)
nr_hvgs <- snakemake@params[["nr_hvgs"]] # how many hvgs for UMAP calc or PCA&UMAP calc
nr_hvgs_batch_correction <- snakemake@params[["nr_hvgs_batch_correction"]] # how many hvgs for batch correction
hvgs_for_batch_correction <- snakemake@params[["hvgs_for_batch_correction"]] # use hvgs?

# rename the already existent Reduced Dimensions
reducedDimNames(sce)[grep("PCA", reducedDimNames(sce))] <- "PCA_before"
reducedDimNames(sce)[grep("UMAP", reducedDimNames(sce))] <- "UMAP_before"

#-------------------------------------------------------------------------------

# calculate hvgs if needed, otherwise keep NULL
if(hvgs_for_batch_correction == TRUE){
  hvgs_for_batch_correction <- modelGeneVar(sce)
  hvgs_for_batch_correction <- getTopHVGs(hvgs_for_batch_correction, 
                                          n=nr_hvgs_batch_correction)
  print(hvgs_for_batch_correction[1:5])
}else{
  hvgs_for_batch_correction <- NULL  
  print("no hvgs")
}

# determine which assay to use based on Rescaling
if(sce$Rescaled[1] == TRUE){
  assay_use <- "corrected"
}else{
  assay_use <- "logcounts"
}

#-------------------------------------------------------------------------------

print("starting MNN")
# run MNNfunction depending on specified method
mnn_fast <- snakemake@params[["mnn_fast"]]
if(mnn_fast == TRUE){
  sce_bc <- correctExperiments(sce, 
                               batch = colData(sce)[,batch_pos], 
                               PARAM = FastMnnParam(),
                               subset.row = hvgs_for_batch_correction, 
                               assay.type = assay_use)
  sce_bc$Correction_method <- rep("FastMNN", ncol(sce_bc))
  
  # calculate UMAP
  hvgs <- modelGeneVar(sce_bc)
  hvgs <- getTopHVGs(hvgs, n=nr_hvgs)
  # can call it PCA because orig PCA is called PCA_before
  reducedDimNames(sce_bc)[reducedDimNames(sce_bc) == "corrected"] <- "PCA"
  colnames(reducedDim(sce_bc, type = "PCA")) <- c(1:ncol(
    reducedDim(sce_bc, type = "PCA")))
  colnames(reducedDim(sce_bc, type = "PCA")) <- paste0(
    "PC", colnames(reducedDim(sce_bc, type = "PCA")))
  
  seeds_umap <- snakemake@params[["seeds_umap"]]
  if(sce_bc$Fraction_ID == "hsc"){
    seed <- seeds_umap[["hsc"]]
  }else if(sce_bc$Fraction_ID == "str"){
    seed <- seeds_umap[["str"]]
  }
  set.seed(seed)
  print(seed)
  sce_bc <- runUMAP(sce_bc, dimred = 'PCA',
                    external_neighbors=TRUE, subset_row = hvgs)
  
}else if(mnn_fast == FALSE){
  sce_bc <- correctExperiments(sce, 
                               batch = colData(sce)[,batch_pos],
                               PARAM = ClassicMnnParam(), 
                               subset.row = hvgs_for_batch_correction,
                               assay.type = assay_use)
  sce_bc$Correction_method <- rep("Classic_MNN", ncol(sce_bc))

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

print(sce_bc)
saveRDS(sce_bc, file = snakemake@output[["sce_09"]])