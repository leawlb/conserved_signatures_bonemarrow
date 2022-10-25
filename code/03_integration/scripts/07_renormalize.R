#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(batchelor, quietly = TRUE)
library(scran)

sce <- readRDS(file = snakemake@input[["sce_06"]])

batch_use <- snakemake@params[["batch_use"]]
rescale <- snakemake@params[["rescale"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# take SCE apart into individuals for quick renormalization
batch_pos <- which(colnames(colData(sce)) == batch_use)
sce$Batch_type_used <- rep(batch_use, ncol(sce))

batch_types <- unique(colData(sce)[,batch_pos])
individual_sce_list <-list()
for(i in batch_types){
  individual_sce_list[[i]] <- sce[,which(colData(sce)[,batch_pos] == i)]
}
print(names(individual_sce_list))

# re-normalize
renorm_list <- multiBatchNorm(individual_sce_list)

# rescale or put back together
print(rescale)
if(rescale == TRUE){
  prepd_sce <- correctExperiments(renorm_list, PARAM = RescaleParam(), assay.type = "logcounts")
  prepd_sce$Rescaled <- rep("TRUE", ncol(prepd_sce))
}else if(rescale == FALSE){
  # merge again
  for(i in 1:length(renorm_list)){
    if(i == 1){ 
      prepd_sce <- renorm_list[[i]]
    }else{
      prepd_sce <- cbind(prepd_sce, renorm_list[[i]])
    }
  }
  prepd_sce$Rescaled <- rep("FALSE", ncol(prepd_sce))
}
print(prepd_sce)

# get hvgs from logcounts only (rescale returns counts in "corrected")
hvgs <- modelGeneVar(prepd_sce)
hvgs <- getTopHVGs(hvgs, n=nr_hvgs)

saveRDS(prepd_sce, file = snakemake@output[["sce_07"]])
saveRDS(hvgs, file = snakemake@output[["hvgs"]])
