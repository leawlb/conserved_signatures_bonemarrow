library(SingleCellExperiment)
library(scuttle)

set.seed(37)

#-------------------------------------------------------------------------------

sce_input_path <- snakemake@input[["sce_input_path"]]
sce_output_path <- snakemake@output[["sce_output_path"]]
print(sce_input_path)
print(sce_output_path)

sce_list <- list()

for(i in 1:length(sce_input_path)){
  sce_list[[i]] <- readRDS(file = sce_input_path[[i]])
}

print(sce_list)

down_output_list <- downsampleBatches(sce_list, assay.type = "counts")

for(i in 1:length(sce_list)){
  
  assays(sce_list[[i]])$downsampled <- down_output_list[[i]]

  stopifnot(c(grepl(sce_list[[i]]$Age_ID[1], sce_output_path[[i]]),
              grepl(sce_list[[i]]$Species_ID[1], sce_output_path[[i]])))
  
  saveRDS(sce_list[[i]], file = sce_output_path[[i]])
}