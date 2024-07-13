#-------------------------------------------------------------------------------
# In this script, transformed DESeq2 objects are generated.
# Hidden sources of variation are taken into account for transformation
# after visual examination which SV to remove from data = add to design

library(DESeq2, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# load object

dsq_list <- base::readRDS(file = snakemake@input[["deseq_input"]])
sva_list <- base::readRDS(file = snakemake@input[["sva"]])

sv_path <- snakemake@input[["sv_path"]]

# info on which SVs to remove from data = put in design
sv_table <- utils::read.csv(
  file = sv_path,
  header = TRUE, 
  sep = ";", 
  check.names=FALSE, 
  stringsAsFactors=FALSE, 
  as.is=TRUE, 
  colClasses = "character")

print(sv_table)

#-------------------------------------------------------------------------------
# filtering out almost empty rows
dsq_list <- lapply(dsq_list, function(dsq){
  
  keep <- BiocGenerics::rowSums(counts(dsq)) >= 10
  dsq <- dsq[keep,]
  return(dsq)
  
})

print(names(dsq_list))

#-------------------------------------------------------------------------------
# add new design with correct number of SVs obtained from SVA (num.sv package)

for(ct in names(dsq_list)){

  print(ct)
  sva <- sva_list[[ct]][[2]]

  # when no SV is found, sva$sv is then 0 = numeric:
  if(class(sva$sv)[1] == "numeric"){
    print("No SVs detected")
    
  # if there are SVs detected, but they shouldn't be removed:
  }else if(!ct %in% sv_table$cell_type){
    print("No SVs to be removed after examination")
    
  # if there are SVs detected and should be removed:
  }else if(ncol(sva$sv) >= 1){
    
    # after examination, which SVs to remove from data -> add to design
    sv_table_temp <- sv_table[sv_table$cell_type == ct,]
    svs_to_remove <- sv_table_temp$SV_to_remove
    svs_to_remove <- base::paste0("sv", svs_to_remove)
    print(base::paste("svs_to_remove:", svs_to_remove))
    
    # prepare sv_df
    sv_df <- base::as.data.frame(sva$sv)
    colnames(sv_df) <- base::paste0("sv", c(1:ncol(sv_df)))
  
    # keep only columns that should be removed = part of the design
    sv_df_temp <- base::as.data.frame(
      sv_df[,c(which(colnames(sv_df) %in% svs_to_remove))])
    # rename to start again with SV1 for compatibility with downstream code
    # !! renamed SVs and svs do NOT necessarily correspond !!
    colnames(sv_df_temp) <- base::paste0("SV", c(1:ncol(sv_df_temp)))
    print(base::paste(
      "remember that SVs and svs do not need to perfectly", 
      "correspond, but the total number must be the same"))
    
    print(head(sv_df))
    print(head(sv_df_temp))
    
    # add SVs to dsq object colData so it can be removed
    colData(dsq_list[[ct]]) <- BiocGenerics::cbind(colData(dsq_list[[ct]]), 
                                                   sv_df_temp)
    
    print(head(colData(dsq_list[[ct]])))
    
    # now, sv_df only contains SVs that should be added to the design
    # this means these SVs will be "removed" from the data
    if(ncol(sv_df_temp) == 1){
      DESeq2::design(dsq_list[[ct]]) <- ~ SV1 + batch + age + condition
    }else if(ncol(sv_df_temp)== 2){
      DESeq2::design(dsq_list[[ct]]) <- ~ SV1 + SV2 + batch + age + condition
    }else if(ncol(sv_df_temp) >= 3){
      stop("number of SVs larger than anticipated")
    }

    print(DESeq2::design(dsq_list[[ct]]))
  }
}

# transform/calculate DGE (normalisation and statistics)
tdsq_list <- lapply(dsq_list, DESeq2::DESeq)

#-------------------------------------------------------------------------------

base::saveRDS(tdsq_list, snakemake@output[["deseq_output"]])

utils::sessionInfo()
