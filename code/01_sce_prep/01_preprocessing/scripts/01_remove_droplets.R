#-------------------------------------------------------------------------------
# author: Amy Danson, Lea Wölbert
# date: 2022-09-20
# remove empty droplets and add doublet score for subsequent doublet removal

library(SingleCellExperiment, quietly = TRUE) 
library(DropletUtils, quietly = TRUE) 
library(scDblFinder, quietly = TRUE) 
set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_input"]])
cutoff_umis <- snakemake@params[["cutoff_umis"]]
cutoff_doublets <-  snakemake@params[["cutoff_doublets"]] 

#-------------------------------------------------------------------------------

# rename rowData, remove duplicated genes
rownames(sce) <- rowData(sce)$Symbol
sce <- sce[!duplicated(rownames(sce)),]

# remove empty droplets
set.seed(37)
out  <- emptyDrops(counts(sce), lower = cutoff_umis)
sce <- sce[,which(out$FDR <= 0.001)]

# identify and remove possible doublets 
set.seed(37)
colData(sce)$doublet_score <- computeDoubletDensity(sce)
sce <- sce[,which(sce$doublet_score <= cutoff_doublets)]

# save
saveRDS(sce, file = snakemake@output[["sce_01"]])