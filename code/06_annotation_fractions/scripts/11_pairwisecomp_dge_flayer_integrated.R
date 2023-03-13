
library(DESeq2)

# in this script, results from DESeq objects are exported for every 
# combination of binary species comparisons
#-------------------------------------------------------------------------------

tdsq_list <- readRDS(file = snakemake@input[["tdsq"]])

a <- "mmus"
b <- "mcas"
c <- "mspr"
d <- "mcar"

comb_list<- list(
  c(a, b),
  c(a, c),
  c(a, d),
  c(b, c),
  c(b, d),
  c(c, d)
)


res_list_list <- lapply(tdsq_list, function(tdsq){
  
  res_list <- lapply(comb_list, function(comb){
 
    print(comb)
    print(comb[1])
    print(comb[2])
  
    res <- results(tdsq, contrast=c("condition", comb[1], comb[2]))
    print(resultsNames(tdsq))
    res_lfcs <- lfcShrink(tdsq, contrast=c("condition", comb[1], comb[2]), type = "ashr")

    return_list <- list(res, res_lfcs)
    names(return_list) <- c("Results", "Results_shrink")
    return(return_list)
  })
  
  for(i in 1:length(comb_list)){
    names(res_list)[[i]] <- paste(comb_list[[i]][1], comb_list[[i]][2], sep = "-")
  }
  print(res_list)
})

names(res_list_list) <- c(1:length(res_list_list))
print(res_list_list)

saveRDS(res_list_list, file = snakemake@output[["res"]])