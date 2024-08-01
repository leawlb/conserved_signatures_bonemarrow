#-------------------------------------------------------------------------------

# in this script, results from DESeq objects are exported for every 
# combination of binary cross-age comparison within species
# for differentially expressed genes (only one combination)

library(DESeq2)

#-------------------------------------------------------------------------------

tdsq_list <- base::readRDS(file = snakemake@input[["deseq_input"]])

fraction_curr <- snakemake@wildcards[["fraction"]]
species <- snakemake@params[["species"]]

celltypes <- names(tdsq_list)
print(celltypes)

#-------------------------------------------------------------------------------
# RESULTS
# get results for young vs. old

res_list <- lapply(tdsq_list, function(tdsq){
  
  res_lfcs <- DESeq2::results(tdsq, contrast=c("age", "yng", "old"))
  
  return(res_lfcs)
})
names(res_list) <- celltypes
names(res_list[[1]])

#-------------------------------------------------------------------------------

base::saveRDS(res_list, file = snakemake@output[["res_list"]])

utils::sessionInfo()
