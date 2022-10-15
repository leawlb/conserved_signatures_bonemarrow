#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(batchelor, quietly = TRUE)

sce <- readRDS(file = snakemake@input[["sce_06"]])

batch_use <- snakemake@params[["batch_use"]]
rescale <- snakemake@params[["rescale"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# take SCE apart into individuals for quick renormalization
individual_sce_list <-list()
for(i in sce[colnames(sce) == batch_use,]){
  sce$Batch_type_used <- rep("batch_use", nrow(sce))
  individual_sce_list[[i]] <- sce[sce$individual == i,]
}

# re-normalize
renorm_list <- multiBatchNorm(individual_sce_list)
# rescale

if(RESCALE == TRUE){
  prepd_sce <- rescaleBatches(renorm_list)
  prepd_sce$Rescaled <- rep("TRUE", nrow(prepd_sce))
}else if(RESCALE == FALSE){
  # merge again
  for(i in 1:length(renorm_list)){
    if(i == 1){ 
      prepd_sce <- renorm_list[[i]]
    }else{
      prepd_sce <- cbind(prepd_sce, renorm_list[[i]])
    }
  }
  prepd_sce$Rescaled <- rep("FALSE", nrow(prepd_sce))
}


# get hvgs 
hvgs <- modelGeneVar(prepd_sce)
hvgs <- getTopHVGs(hvgs, n=2000)
