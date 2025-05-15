#-------------------------------------------------------------------------------
# separate SCE based on cell fate probabilities
# the cut-off for including cells into a branch is calculated from
# the cell fate probability distribution in "Early MPPs", as the "last" cell 
# type from the early compartment
# The cut-off will be set to the first quartile of cell fate probabilities for
# each lineage -1.5 * IQR 
# this will only remove few outliers from Early MPPs from each branch (and)
# cell with values below the cut-off

set.seed(37)
library(SingleCellExperiment)

#-------------------------------------------------------------------------------
# load data
sce_hsc <- base::readRDS(snakemake@input[["sce_hsc"]])
sce_pseudotime <- base::readRDS(snakemake@input[["sce_pseudotime"]])

#-------------------------------------------------------------------------------
# transfer the data from the pseudotime object to the sce object with
# correct UMAP coordinates and other metadata
sce_hsc$pseudotime <- sce_pseudotime$dpt_pseudotime

# know the correct sequence from the adata object
sce_hsc$Lymphoid <-  SingleCellExperiment::reducedDims(sce_pseudotime)$lineages_fwd[,1]
sce_hsc$Erythroid <-  SingleCellExperiment::reducedDims(sce_pseudotime)$lineages_fwd[,2]
sce_hsc$Neutrophil <-  SingleCellExperiment::reducedDims(sce_pseudotime)$lineages_fwd[,3]

#-------------------------------------------------------------------------------
# get all cell fate probability values for each lineage from cell type
vals_ery <- sce_hsc$Erythroid[sce_hsc$celltypes == "Early MPP"]
# calculate the cutoff as described
cut_off_ery <- base::summary(vals_ery)["1st Qu."] - 1.5*stats::IQR(vals_ery)

vals_lym <- sce_hsc$Lymphoid[sce_hsc$celltypes == "Early MPP"]
cut_off_lym <- base::summary(vals_lym)["1st Qu."] - 1.5*stats::IQR(vals_lym)

vals_neu <- sce_hsc$Neutrophil[sce_hsc$celltypes == "Early MPP"]
cut_off_neu <- base::summary(vals_neu)["1st Qu."] - 1.5*stats::IQR(vals_neu)

# separate
sce_ery <- sce_hsc[,sce_hsc$Erythroid > cut_off_ery]
sce_lym <- sce_hsc[,sce_hsc$Lymphoid > cut_off_lym]
sce_neu <- sce_hsc[,sce_hsc$Neutrophil > cut_off_neu]

#-------------------------------------------------------------------------------
# save data
base::saveRDS(sce_ery, snakemake@output[["sce_ery"]])
base::saveRDS(sce_lym, snakemake@output[["sce_lym"]])
base::saveRDS(sce_neu, snakemake@output[["sce_neu"]])

utils::sessionInfo()
