
library(DESeq2)
set.seed(37)

# In this script, Transformed DESEQ objects with two different designs are
# generated. One design takes SVs into account (After), the other doesn't (BEFORE)

#-------------------------------------------------------------------------------
# load object

dsq_list <- readRDS(file = snakemake@input[["dsq_11"]])
sva_list <- readRDS(file = snakemake@input[["sva_11"]])

#-------------------------------------------------------------------------------
# filtering out almost empty rows
dsq_list <- lapply(dsq_list, function(dsq){
  
  keep <- rowSums(counts(dsq)) >= 10
  dsq <- dsq[keep,]
  return(dsq)
})
print(dsq_list)
#-------------------------------------------------------------------------------
# transform
tdsq_list <- lapply(dsq_list, DESeq)
saveRDS(tdsq_list, snakemake@output[["tdsq_before_11"]])

#-------------------------------------------------------------------------------
# add new design with correct number of SVs obtained from SVA
for(i in 1:length(dsq_list)){

  sva <- sva_list[[i]][[2]]
  print(sva$sv)
  
  if(sva$sv == 0){
    print("No hidden variation: no design change required")
  }else if(ncol(sva$sv) >= 1){
    for(j in 1:ncol(sva$sv)){
      colData(dsq_list[[i]])[,ncol(colData(dsq_list[[i]]))+1] <- sva$sv[,j]
      colnames(colData(dsq_list[[i]]))[ncol(colData(dsq_list[[i]]))] <- paste0("SV", j)
    }
  
    if(ncol(sva$sv) == 1){
      design(dsq_list[[i]]) <- ~ SV1 + batch + age + condition
    }else if(ncol(sva$sv) == 2){
      design(dsq_list[[i]]) <- ~ SV1 + SV2 + batch + age + condition
    }else if(ncol(sva$sv) == 3){
      design(dsq_list[[i]]) <- ~ SV1 + SV2 + SV3 + batch + age + condition
    }else if(ncol(sva$sv) == 4){
      design(dsq_list[[i]]) <- ~ SV1 + SV2 + SV3 + SV4 + batch + age + condition
    }else if(ncol(sva$sv) >= 5){
      stop("number of SVs larger than anticipated")
    }
    print(colData(dsq_list[[i]]))
    print(i)
    print(design(dsq_list[[i]]))
  }
}

# rerun DEseq with new design and coldata
tdsq_list_after <- lapply(dsq_list, DESeq)
saveRDS(tdsq_list_after, snakemake@output[["tdsq_after_11"]])
