#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(batchelor, quietly = TRUE)

# load object and prepare batch vector for correctExperiments()
print(snakemake@input[["sce_08"]])
sce <- readRDS(file = snakemake@input[[1]])
batch_use <- snakemake@params[["batch_use"]]
batch_pos <- which(colnames(colData(sce)) == batch_use)
nr_hvgs <- snakemake@params[["nr_hvgs"]]
hvgs_for_batch_correction <- snakemake@params[["hvgs_for_batch_correction"]]

source(file = snakemake@params[["sce_functions"]])

# rename the already existent Reduced Dimensions
reducedDimNames(sce)[grep("PCA", reducedDimNames(sce))] <- "PCA_before"
reducedDimNames(sce)[grep("UMAP", reducedDimNames(sce))] <- "UMAP_before"

# load hvgs if needed, otherwise keep NULL
if(hvgs_for_batch_correction == TRUE){
  hvgs <- readRDS(file = snakemake@input[[2]])
  print(hvgs[1:10])
}else{
  hvgs <- NULL  
}

# determine which assay to use based on Rescaling
if(sce$Rescaled[1] == TRUE){
  assay_use <- "corrected"
}else{
  assay_use <- "logcounts"
}

print("starting MNN")
# run MNNfunction depending on specified method
mnn_fast <- snakemake@params[["mnn_fast"]]
if(mnn_fast == TRUE){
  sce_bc <- correctExperiments(sce, 
                               batch = colData(sce)[,batch_pos], 
                               PARAM = FastMnnParam(),
                               subset.row = hvgs, 
                               assay.type = assay_use)
  sce_bc$Correction_method <- rep("FastMNN", ncol(sce_bc))
}else if(mnn_fast == FALSE){
  sce_bc <- correctExperiments(sce, 
                               batch = colData(sce)[,batch_pos],
                               PARAM = ClassicMnnParam(), 
                               subset.row = hvgs,
                               assay.type = assay_use)
  sce_bc$Correction_method <- rep("Classic_MNN", ncol(sce_bc))
}

# Reduce dimensions
sce <- reduce_dims(sce, nr_hvgs = nr_hvgs)

print(sce_bc)

saveRDS(sce_bc, file = snakemake@output[["sce_09"]])
