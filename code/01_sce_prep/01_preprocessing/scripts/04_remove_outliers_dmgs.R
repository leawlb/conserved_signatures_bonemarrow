#-------------------------------------------------------------------------------
# authors: Amy Danson, Lea Wölbert
# remove cells with low quality, remove dmgs

library(scuttle, quietly = TRUE) 
set.seed(37)

sce <- base::readRDS(file = snakemake@input[["sce_input"]])
dmg_list <- base::readRDS(file = snakemake@input[["dmg_list"]])

# QC cutoffs 
cutoff_sum <- snakemake@params[["cutoff_sum"]] 
cutoff_detected <- snakemake@params[["cutoff_detected"]]
cutoff_mitos <- snakemake@params[["cutoff_mitos"]]

#-------------------------------------------------------------------------------

# remove all dmgs
print(sce)
sce <- sce[!rownames(sce) %in% dmg_list,]
print(sce)

#-------------------------------------------------------------------------------

# get percentage of mitochondrial genes
mito_genes <- base::grep("mt-", rownames(sce))
print(mito_genes)
qcdf<- scuttle::perCellQCMetrics(sce, subsets=list(Mito=mito_genes)) 

# identify and remove cells below QC cutoffs 
sum_out <- qcdf$sum < cutoff_sum 
det_out <- qcdf$detected < cutoff_detected 
mito_out <- qcdf$subsets_Mito_percent > cutoff_mitos

remove_pos <- sum_out | det_out | mito_out

print(sce)
sce <- sce[,!remove_pos]
print(sce)

base::saveRDS(sce, file = snakemake@output[["sce_output"]])

utils::sessionInfo()