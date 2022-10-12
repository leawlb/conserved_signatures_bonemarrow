#-------------------------------------------------------------------------------
# author: Amy Danson, Lea Wölbert
# date: 2022-09-20
# remove outlier cells with low quality (low reads, small library, high mito)
# evaluate in report.Rmd

library(SingleCellExperiment, quietly = TRUE) 
library(scuttle, quietly = TRUE) 

sce <- readRDS(file = snakemake@input[["sce_02"]])

# calculate quality control dataframe
qcdf <- perCellQCMetrics(sce)

# set QC cutoffs below which cells are removed
cutoff_sum <- snakemake@params[["cutoff_sum"]] 
cutoff_detected <- snakemake@params[["cutoff_detected"]]

sum_out <- qcdf$sum < cutoff_sum 
det_out <- qcdf$sum < cutoff_detected 

# remove all cells to discard and add info to qc
remove_pos <- sum_out | det_out 
sce <- sce[,!remove_pos]

# save as new objects
saveRDS(sce, file = snakemake@output[["sce_03"]])