#-------------------------------------------------------------------------------
# in this script, results from DESeq objects are exported for every 
# combination of binary species comparisons
# for differentially expressed genes

library(DESeq2)
library(tidyverse)
#library(ashr)

#-------------------------------------------------------------------------------

tdsq_list <- readRDS(file = snakemake@input[["deseq_input"]])

fraction_curr <- snakemake@wildcards[["fraction"]]
species <- snakemake@params[["species"]]

celltypes <- names(tdsq_list)
print(celltypes)

#-------------------------------------------------------------------------------

a <- species[1]
b <- species[2]
c <- species[3]
d <- species[4]

# every combination for easier sub-setting later
comb_list<- list(
  c(a, b),
  c(a, c),
  c(a, d),
  c(b, a),
  c(b, c),
  c(b, d),
  c(c, a),
  c(c, b),
  c(c, d),
  c(d, a),
  c(d, b),
  c(d, c)
)

#-------------------------------------------------------------------------------
# RESULTS
# get results for each combination

res_list_list <- lapply(tdsq_list, function(tdsq){
  
  res_list <- lapply(comb_list, function(comb){
    
    res_lfcs <- DESeq2::results(tdsq, contrast=c("species", comb[1], comb[2]))
    
    #print(res_lfcs)
    return(res_lfcs)
  })
  for(i in 1:length(comb_list)){
    names(res_list)[[i]] <- paste0(comb_list[[i]][1], "-", comb_list[[i]][2])
  }
  return(res_list)
})
names(res_list_list) <- celltypes
names(res_list_list[[1]])

saveRDS(res_list_list, file = snakemake@output[["celltype_res"]])

#-------------------------------------------------------------------------------
# Get results DF

# add info
res_df_list <- lapply(res_list_list, function(res_list){
  
  res_df <- data.frame()
  for(j in names(res_list)){
    res_temp <- data.frame(res_list[[j]])
    res_temp <- rownames_to_column(res_temp, var = "gene")

    res_temp$comparison <- j
    res_temp$species <- str_sub(j, 1, 4)

    res_df <- rbind(res_df, res_temp)
  }
  return(res_df)
})
names(res_df_list) <- names(res_list_list)

#-------------------------------------------------------------------------------
saveRDS(res_df_list, file = snakemake@output[["celltype_res_dfs"]])

sessionInfo()
