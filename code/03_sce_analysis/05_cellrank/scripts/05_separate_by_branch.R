#-------------------------------------------------------------------------------
# separate SCE based on cell fate probabilities
# the cut-off for including cells into a branch is calculated from
# the cell fate probability distribution in "Early MPPs", as the "last" cell 
# type from the early compartment
# The cut-off will be set to the first quartile of cell fate probabilities for
# each lineage -1.5 * IQR 
# this will only remove few outliers from HSCs and MPPs from each branch

set.seed(37)

#-------------------------------------------------------------------------------
# load data
sce_hsc <- readRDS(snakemake@input[["sce_hsc"]])
sce_pseudotime <- readRDS(snakemake@input[["sce_pseudotime"]])

#-------------------------------------------------------------------------------
# transfer the data from the pseudotime object to the sce object with
# correct UMAP coordinates and other metadata
sce_hsc$pseudotime <- sce_pseudotime$dpt_pseudotime

# know the correct sequence from the adata object
sce_hsc$Erythroid <-  reducedDims(sce_pseudotime)$lineages_fwd[,1]
sce_hsc$Lymphoid <-  reducedDims(sce_pseudotime)$lineages_fwd[,2]
sce_hsc$Neutrophil <-  reducedDims(sce_pseudotime)$lineages_fwd[,3]

#-------------------------------------------------------------------------------
# get all cell fate probability values for each lineage from cell type
vals_ery <- sce_hsc$Erythroid[sce_hsc$celltypes == "Early MPPs"]
# calculate the cutoff as described
cut_off_ery <- summary(vals_ery)["1st Qu."] - 1.5*IQR(vals_ery)

vals_lym <- sce_hsc$Lymphoid[sce_hsc$celltypes == "Early MPPs"]
cut_off_lym <- summary(vals_lym)["1st Qu."] - 1.5*IQR(vals_lym)

vals_neu <- sce_hsc$Neutrophil[sce_hsc$celltypes == "Early MPPs"]
cut_off_neu <- summary(vals_neu)["1st Qu."] - 1.5*IQR(vals_neu)

# separate
sce_ery <- sce_hsc[,sce_hsc$Erythroid > cut_off_ery]
sce_lym <- sce_hsc[,sce_hsc$Lymphoid > cut_off_lym]
sce_neu <- sce_hsc[,sce_hsc$Neutrophil > cut_off_neu]

#-------------------------------------------------------------------------------
# save data
saveRDS(sce_ery, snakemake@output[["sce_ery"]])
saveRDS(sce_lym, snakemake@output[["sce_lym"]])
saveRDS(sce_neu, snakemake@output[["sce_neu"]])
