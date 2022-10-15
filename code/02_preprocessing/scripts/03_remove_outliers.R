#-------------------------------------------------------------------------------
# author: Amy Danson, Lea Wölbert
# date: 2022-09-20
# remove outlier cells with low quality (low reads, small library, high mito)
# evaluate in report.Rmd

library(SingleCellExperiment, quietly = TRUE) 
library(scuttle, quietly = TRUE) 

sce <- readRDS(file = snakemake@input[["sce_02"]])

# set QC cutoffs below which cells are removed
cutoff_sum <- snakemake@params[["cutoff_sum"]] 
cutoff_detected <- snakemake@params[["cutoff_detected"]]
cutoff_mitos <- snakemake@params[["cutoff_mitos"]]

# get mitos
mito_genes <- grep("mt-", rownames(sce))
qcdf<- perCellQCMetrics(sce, 
                        subsets=list(Mito=mito_genes)) # before outlier removal

# set QC cutoffs below which cells are removed
sum_out <- qcdf$sum < cutoff_sum 
det_out <- qcdf$detected < cutoff_detected 
mito_out <- qcdf$subsets_Mito_percent > cutoff_mitos

# remove all cells to discard and add info to qc
remove_pos <- sum_out | det_out | mito_out
qcdf$removed <- remove_pos

sce <- sce[,!remove_pos]

# save as new objects
saveRDS(sce, file = snakemake@output[["sce_03"]])