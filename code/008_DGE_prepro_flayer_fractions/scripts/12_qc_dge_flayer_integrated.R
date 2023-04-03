
library(DESeq2)
library(sva)
set.seed(37)

# In this script, QC data is exported, including logr-transformed objects
# for visualisation only, and detected hidden sources of variation (SVs)

#-------------------------------------------------------------------------------
# load object

dsq_list <- readRDS(file = snakemake@input[["dsq_11"]])

#-------------------------------------------------------------------------------
# filter out almost empty rows
dsq_list <- lapply(dsq_list, function(dsq){
  
  keep <- rowSums(counts(dsq)) >= 10
  dsq <- dsq[keep,]
  return(dsq)
})
print(dsq_list)
#-------------------------------------------------------------------------------
# transform counts FOR VISUALISATION ONLY
# log transformation and sequencing depth correction 
rld_list <- lapply(dsq_list, rlog, blind = FALSE)
print(rld_list)
saveRDS(rld_list, snakemake@output[["rld_11"]])

#-------------------------------------------------------------------------------
# get the number of hidden sources of variations

# transformation is required for num.sv
tdsq_list <- lapply(dsq_list, DESeq)

sva_list <- lapply(tdsq_list, function(dsq){
  
  data <- counts(dsq, normalized = TRUE)
  data  <- data[which(rowMeans(data) > 2), ]
  
  mod <- model.matrix(~ batch + age + condition, colData(dsq))
  mod0 <- model.matrix(~ 1, colData(dsq))
  
  n_sv_be <- num.sv(data, mod, method = "be")  
  n_sv_lk <- num.sv(data, mod, method = "leek")  
  
  # n_sv_lk is so large that it usually doesn't work 
  
  svseq <- svaseq(data, mod, mod0, n.sv = n_sv_be)
  print(c(n_sv_be, n_sv_lk))
  return(list(c(n_sv_be, n_sv_lk), svseq))
})

saveRDS(sva_list, snakemake@output[["sva_11"]])
