#-------------------------------------------------------------------------------

library(monocle3)

set.seed(37)

#-------------------------------------------------------------------------------

cds <- readRDS(snakemake@input[["cds_input"]])

# trajectory
cds <- learn_graph(cds, use_partition = TRUE)

# save image to select root nodes
#pdf(snakemake@output[["pdf_output"]])
#plot_cells_3d(cds, label_principal_points = TRUE, color_cells_by = "partition")
#dev.off() 

saveRDS(cds, snakemake@output[["cds_output"]])
