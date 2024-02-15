#-------------------------------------------------------------------------------
library(SingleCellExperiment)
set.seed(37)

sce_path_inp <- snakemake@input[["sce_input"]]
sce_path_out <- snakemake@output[["sce_output"]]

#-------------------------------------------------------------------------------
# read objects

sce_hsc <- readRDS(file = sce_path_inp[[1]])
sce_str <- readRDS(file = sce_path_inp[[2]])

#-------------------------------------------------------------------------------

# rename celltypes to remove spaces
levels_hsc_new <- gsub(" ", "_", levels(sce_hsc$celltypes))
levels_hsc_new <- gsub("/", "_", levels_hsc_new)
levels_str_new <- gsub(" ", "_", levels(sce_str$celltypes))
levels_str_new <- gsub("/", "_", levels_str_new)
sce_hsc$celltypes_copy <- sce_hsc$celltypes
sce_str$celltypes_copy <- sce_str$celltypes
sce_hsc$celltypes <- gsub(" ", "_", sce_hsc$celltypes)
sce_hsc$celltypes <- gsub("/", "_", sce_hsc$celltypes)
sce_str$celltypes <- gsub(" ", "_", sce_str$celltypes)
sce_str$celltypes <- gsub("/", "_", sce_str$celltypes)
sce_hsc$celltypes <- factor(sce_hsc$celltypes, levels_hsc_new)
sce_str$celltypes <- factor(sce_str$celltypes, levels_str_new)
# separate objects

print(levels(sce_hsc$celltypes))
print(levels(sce_str$celltypes))

sce_hsc_list <- list()
for(i in levels(sce_hsc$celltypes)){
  sce_hsc_list[[i]] <- sce_hsc[,sce_hsc$celltypes == i]
}

sce_str_list <- list()
for(i in levels(sce_str$celltypes)){
  sce_str_list[[i]] <- sce_str[,sce_str$celltypes == i]
}

#-------------------------------------------------------------------------------
# save separated objects

sce_path_out_hsc <- sce_path_out[grep("hsc", sce_path_out)]
sce_path_out_str <- sce_path_out[grep("str", sce_path_out)]


print(sce_path_out_hsc)
print(levels(sce_hsc$celltypes))
for(i in levels(sce_hsc$celltypes)){
  print(i)
  output_path <- sce_path_out_hsc[[grep(i, sce_path_out_hsc)]]
  print(output_path)
  saveRDS(sce_hsc_list[[i]], file = output_path)
}
print("finished HSCs")
print(sce_path_out_str) 
for(i in levels(sce_str$celltypes)){
  print(i)
  output_path <- sce_path_out_str[[grep(i, sce_path_out_str)]]
  print(output_path)
  saveRDS(sce_str_list[[i]], file = output_path)
}
