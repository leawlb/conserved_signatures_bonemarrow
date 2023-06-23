#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(batchelor, quietly = TRUE)
library(scran, quietly = TRUE)
set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_input"]])

batch_use <- snakemake@params[["batch_use"]]
rescale <- snakemake@params[["rescale"]]

# take SCE apart, into individuals for quick renormalization
batch_pos <- which(colnames(colData(sce)) == batch_use)
sce$Batch_type_used <- rep(batch_use, ncol(sce))

batch_types <- unique(colData(sce)[,batch_pos])
individual_sce_list <-list()
for(i in batch_types){
  individual_sce_list[[i]] <- sce[,which(colData(sce)[,batch_pos] == i)]
}
print(names(individual_sce_list))

individual_sce_list <- lapply(individual_sce_list, function(sce){
  assays(sce)$logcounts_orig <- assays(sce)$logcounts
  return(sce)
})

# re-normalize
renorm_list <- multiBatchNorm(individual_sce_list)

# rescale or put back together
print(rescale)
if(rescale == TRUE){
  prepd_sce <- correctExperiments(renorm_list, PARAM = RescaleParam(), 
                                  assay.type = "logcounts")
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

saveRDS(prepd_sce, file = snakemake@output[["sce_01"]])