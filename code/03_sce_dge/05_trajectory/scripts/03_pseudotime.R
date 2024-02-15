library(monocle3)

set.seed(37)

#-------------------------------------------------------------------------------

# prepare SCE object for conversion
cds <- readRDS(snakemake@input[["cds_input"]])

#nodes <- c("Y_1400", "Y_1274", "Y_1362", "Y_995")
nodes <- c("Y_1400", "Y_1274", "Y_1362", "Y_1777")
cds <- order_cells(cds, root_pr_nodes = nodes)

#plot_cells(cds, color_cells_by = "pseudotime")
#plot_cells(cds, color_cells_by = "celltypes")

saveRDS(cds, snakemake@output[["cds_output"]])
