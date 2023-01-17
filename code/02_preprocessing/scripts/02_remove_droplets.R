#-------------------------------------------------------------------------------
# author: Amy Danson, Lea Wölbert
# date: 2022-09-20
# remove empty droplets and add doublet score for subsequent doublet removal

print(getwd())

library(SingleCellExperiment, quietly = TRUE) 
library(DropletUtils, quietly = TRUE) 
library(scDblFinder, quietly = TRUE) 
set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_01"]])

# change rownames to Symbols for better readability
rownames(sce) <- rowData(sce)$Symbol
# remove duplicated rows
sce <- sce[!duplicated(rownames(sce)),]


# set a cutoff of minimum UMI codes required for a barcode to count as a cell
# evaluate in report.Rmd
cutoff_umis <- snakemake@params[["cutoff_umis"]]

set.seed(37)
# remove barcodes that are likely empty
out  <- emptyDrops(counts(sce), lower = cutoff_umis)
sce <- sce[,which(out$FDR <= 0.001)]

# identify possible doublets and add to coldata
colData(sce)$doublet_score <- computeDoubletDensity(sce)

# set a doublet_score limit above which possible doublets are excluded
# evaluate in report.Rmd
# setting a very high limit so no barcodes are accidentally removed
set.seed(37)
cutoff_doublets <-  snakemake@params[["cutoff_doublets"]] 
sce <- sce[,which(sce$doublet_score <= cutoff_doublets)]

# save as new object
saveRDS(sce, file = snakemake@output[["sce_02"]])