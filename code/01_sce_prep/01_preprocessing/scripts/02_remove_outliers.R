#-------------------------------------------------------------------------------
# author: Amy Danson, Lea Wölbert
# date: 2022-09-20
# remove cells with low quality

library(SingleCellExperiment, quietly = TRUE) 
library(scuttle, quietly = TRUE) 
set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_01"]])

# QC cutoffs 
cutoff_sum <- snakemake@params[["cutoff_sum"]] 
cutoff_detected <- snakemake@params[["cutoff_detected"]]
cutoff_mitos <- snakemake@params[["cutoff_mitos"]]

#-------------------------------------------------------------------------------

# get mitos
mito_genes <- grep("mt-", rownames(sce))
print(mito_genes)
qcdf<- perCellQCMetrics(sce, 
                        subsets=list(Mito=mito_genes)) # before outlier removal

# identify and remove cells below QC cutoffs 
sum_out <- qcdf$sum < cutoff_sum 
det_out <- qcdf$detected < cutoff_detected 
mito_out <- qcdf$subsets_Mito_percent > cutoff_mitos

remove_pos <- sum_out | det_out | mito_out

sce <- sce[,!remove_pos]

# save 
saveRDS(sce, file = snakemake@output[["sce_02"]])