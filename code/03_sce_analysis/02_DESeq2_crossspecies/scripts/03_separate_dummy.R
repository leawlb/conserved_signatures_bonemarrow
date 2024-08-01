#-------------------------------------------------------------------------------
set.seed(37)

sce_path_inp <- snakemake@input[["sce_input"]]
path_out <- snakemake@output[["output"]]

#-------------------------------------------------------------------------------
# read objects for getting the celltypes

sce_hsc <- base::readRDS(file = sce_path_inp[[1]])
sce_str <- base::readRDS(file = sce_path_inp[[2]])

#-------------------------------------------------------------------------------

# rename celltypes to remove spaces
levels_hsc_new <- base::gsub(" ", "_", levels(sce_hsc$celltypes))
levels_hsc_new <- base::gsub("/", "_", levels_hsc_new)
levels_hsc_new <- base::gsub("[.]", "", levels_hsc_new)
levels_str_new <- base::gsub(" ", "_", levels(sce_str$celltypes))
levels_str_new <- base::gsub("/", "_", levels_str_new)
levels_str_new <- base::gsub("[.]", "", levels_str_new)

sce_hsc$celltypes_copy <- sce_hsc$celltypes
sce_str$celltypes_copy <- sce_str$celltypes

sce_hsc$celltypes <- base::gsub(" ", "_", sce_hsc$celltypes)
sce_hsc$celltypes <- base::gsub("/", "_", sce_hsc$celltypes)
sce_hsc$celltypes <- base::gsub("[.]", "", sce_hsc$celltypes)
sce_str$celltypes <- base::gsub(" ", "_", sce_str$celltypes)
sce_str$celltypes <- base::gsub("/", "_", sce_str$celltypes)
sce_str$celltypes <- base::gsub("[.]", "", sce_str$celltypes)

print(levels_hsc_new)
print(levels_str_new)

sce_hsc$celltypes <- factor(sce_hsc$celltypes, levels_hsc_new)
sce_str$celltypes <- factor(sce_str$celltypes, levels_str_new)

#-------------------------------------------------------------------------------
# save just the cell type as

path_out_hsc <- path_out[base::grep("hsc", path_out)]
path_out_str <- path_out[base::grep("str", path_out)]

print(path_out_hsc)
print(levels(sce_hsc$celltypes))
for(ct in levels(sce_hsc$celltypes)){
  print(ct)
  output_path <- path_out_hsc[[base::grep(ct, path_out_hsc)]]
  print(output_path)
  base::saveRDS(ct, file = output_path)
}

print(path_out_str) 
for(ct in levels(sce_str$celltypes)){
  print(ct)
  output_path <- path_out_str[[base::grep(ct, path_out_str)]]
  print(output_path)
  base::saveRDS(ct, file = output_path)
}

utils::sessionInfo()