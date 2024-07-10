#-------------------------------------------------------------------------------
# authors: Amy Danson, Lea Wölbert
# remove empty dropletsm add doublet score for doublet removal, remove doublets
# remove duplicated genes

library(DropletUtils, quietly = TRUE) 
library(scDblFinder, quietly = TRUE) 
library(SingleCellExperiment, quietly = TRUE) 
set.seed(37)

# load objects and params defined in corresponding snakefile
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

cutoff_umis <- snakemake@params[["cutoff_umis"]]
cutoff_doublets <- snakemake@params[["cutoff_doublets"]] 

#-------------------------------------------------------------------------------

# rename rowData, remove duplicated genes
rownames(sce) <- rowData(sce)$Symbol
sce <- sce[!base::duplicated(rownames(sce)),]

# remove empty droplets
out  <- DropletUtils::emptyDrops(counts(sce), 
                                 lower = cutoff_umis)
sce <- sce[,which(out$FDR <= 0.001)]

# identify and remove possible doublets 
colData(sce)$doublet_score <- scDblFinder::computeDoubletDensity(sce)
sce <- sce[,which(sce$doublet_score <= cutoff_doublets)]

# save
base::saveRDS(sce, file = snakemake@output[["sce_output"]])

# print session info
utils::sessionInfo()