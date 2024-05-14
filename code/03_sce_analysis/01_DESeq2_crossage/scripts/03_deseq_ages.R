#-------------------------------------------------------------------------------

# In this script, Transformed DESeq2 objects are generated.
# Hidden sources of variation are taken into account for transformation

library(DESeq2)
set.seed(37)

#-------------------------------------------------------------------------------
# load object

dsq_list <- readRDS(file = snakemake@input[["deseq_input"]])
sva_list <- readRDS(file = snakemake@input[["sva"]])

#-------------------------------------------------------------------------------
# filtering out almost empty rows
dsq_list <- lapply(dsq_list, function(dsq){
  
  keep <- rowSums(counts(dsq)) >= 10
  dsq <- dsq[keep,]
  return(dsq)
})
print(dsq_list)

#-------------------------------------------------------------------------------
# add new design with correct number of SVs obtained from SVA
for(i in 1:length(dsq_list)){

  sva <- sva_list[[i]][[2]]
  print(sva$sv)

  # sometimes, no SV is found, sva$sv is then 0
  if(class(sva$sv)[1] == "numeric"){
    print("No hidden variation: no design change required")
  }else if(ncol(sva$sv) >= 1){
    
    for(j in 1:ncol(sva$sv)){
      colData(dsq_list[[i]])[,ncol(colData(dsq_list[[i]]))+1] <- sva$sv[,j]
      colnames(colData(dsq_list[[i]]))[
        ncol(colData(dsq_list[[i]]))] <- paste0("SV", j)
    }
  
    if(ncol(sva$sv) == 1){
      design(dsq_list[[i]]) <- ~ SV1 + age 
    }else if(ncol(sva$sv) == 2){
      design(dsq_list[[i]]) <- ~ SV1 + SV2 + age 
    }else if(ncol(sva$sv) == 3){
      design(dsq_list[[i]]) <- ~ SV1 + SV2 + SV3 + age 
    }else if(ncol(sva$sv) == 4){
      design(dsq_list[[i]]) <- ~ SV1 + SV2 + SV3 + SV4 + age 
    }else if(ncol(sva$sv) >= 5){
      stop("number of SVs larger than anticipated")
    }
    print(i)
    print(design(dsq_list[[i]]))
  }
}

# transform/calculate DGE (normalisation and statistics)
tdsq_list <- lapply(dsq_list, DESeq2::DESeq)

#-------------------------------------------------------------------------------

saveRDS(tdsq_list, snakemake@output[["deseq_output"]])

sessionInfo()