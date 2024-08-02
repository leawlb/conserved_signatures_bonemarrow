#-------------------------------------------------------------------------------
# prepare for DESeq2 QC 
# In this script, QC data is exported, including logr-transformed objects
# for visualisation only, as well as detected hidden sources of variation (SVs)

library(DESeq2, quietly = TRUE)
library(sva, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# load object

dsq_list <- base::readRDS(file = snakemake@input[["deseq_input"]])

#-------------------------------------------------------------------------------
# filter out almost empty rows
dsq_list <- lapply(dsq_list, function(dsq){
  
  keep <- BiocGenerics::rowSums(BiocGenerics::counts(dsq)) >= 10
  dsq <- dsq[keep,]
  return(dsq)
  
})
print(dsq_list)

#-------------------------------------------------------------------------------
# transform and save counts FOR VISUALISATION ONLY
# log transformation and sequencing depth correction, basically normalization

rld_list <- lapply(dsq_list, DESeq2::rlog, blind = TRUE)
print(rld_list)

base::saveRDS(rld_list, snakemake@output[["rlog"]])

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# get the number of hidden sources of variations
# taken from this tutorial:
# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

# transformation is required for num.sv
tdsq_list <- lapply(dsq_list, DESeq2::DESeq)

# get a list of dataframes containing sv data
sva_list <- lapply(tdsq_list, function(dsq){
  
  data <- BiocGenerics::counts(dsq, normalized = TRUE)
  data <- data[which(BiocGenerics::rowMeans(data) > 1), ] # like tutorial
  
  # full model matrix
  mod <- stats::model.matrix(~ batch + age + species, colData(dsq))
  # null model matrix
  mod0 <- stats::model.matrix(~ 1, colData(dsq))
  
  # try two different methods but use BE
  n_sv_be <- sva::num.sv(data, mod, method = "be")  
  n_sv_lk <- sva::num.sv(data, mod, method = "leek")  

  svseq <- sva::svaseq(data, mod, mod0, n.sv = n_sv_be)
  print(c(n_sv_be, n_sv_lk))
  
  return(list(c(n_sv_be, n_sv_lk), svseq))
})

#-------------------------------------------------------------------------------
base::saveRDS(sva_list, snakemake@output[["sva"]])

utils::sessionInfo()
