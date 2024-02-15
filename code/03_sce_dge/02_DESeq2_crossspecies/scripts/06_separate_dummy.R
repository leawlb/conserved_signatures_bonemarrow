#-------------------------------------------------------------------------------
library(SingleCellExperiment)
set.seed(37)

sce_path_inp <- snakemake@input[["sce_input"]]
path_out <- snakemake@output[["output"]]

#-------------------------------------------------------------------------------
# read objects for getting the celltypes

sce_hsc <- readRDS(file = sce_path_inp[[1]])
sce_str <- readRDS(file = sce_path_inp[[2]])

#-------------------------------------------------------------------------------

# rename celltypes to remove spaces
levels_hsc_new <- gsub(" ", "_", levels(sce_hsc$celltypes))
levels_hsc_new <- gsub("/", "_", levels_hsc_new)
levels_hsc_new <- gsub("[.]", "", levels_hsc_new)
levels_str_new <- gsub(" ", "_", levels(sce_str$celltypes))
levels_str_new <- gsub("/", "_", levels_str_new)
levels_str_new <- gsub("[.]", "", levels_str_new)
sce_hsc$celltypes_copy <- sce_hsc$celltypes
sce_str$celltypes_copy <- sce_str$celltypes
sce_hsc$celltypes <- gsub(" ", "_", sce_hsc$celltypes)
sce_hsc$celltypes <- gsub("/", "_", sce_hsc$celltypes)
sce_hsc$celltypes <- gsub("[.]", "", sce_hsc$celltypes)
sce_str$celltypes <- gsub(" ", "_", sce_str$celltypes)
sce_str$celltypes <- gsub("/", "_", sce_str$celltypes)
sce_str$celltypes <- gsub("[.]", "", sce_str$celltypes)

print(levels_hsc_new)
print(levels_str_new)


sce_hsc$celltypes <- factor(sce_hsc$celltypes, levels_hsc_new)
sce_str$celltypes <- factor(sce_str$celltypes, levels_str_new)


#-------------------------------------------------------------------------------
# save separated objects

path_out_hsc <- path_out[grep("hsc", path_out)]
path_out_str <- path_out[grep("str", path_out)]

print(path_out_hsc)
print(levels(sce_hsc$celltypes))
for(i in levels(sce_hsc$celltypes)){
  print(i)
  output_path <- path_out_hsc[[grep(i, path_out_hsc)]]
  print(output_path)
  saveRDS(i, file = output_path)
}

print(path_out_str) 
for(i in levels(sce_str$celltypes)){
  print(i)
  output_path <- path_out_str[[grep(i, path_out_str)]]
  print(output_path)
  saveRDS(i, file = output_path)
}
