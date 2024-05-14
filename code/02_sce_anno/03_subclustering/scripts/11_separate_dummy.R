#-------------------------------------------------------------------------------

set.seed(37)

#-------------------------------------------------------------------------------
# read objects

sce_path_inp <- snakemake@input[["sce_input"]]
path_out <- snakemake@output[["output"]]

sce_hsc <- readRDS(file = sce_path_inp[[1]])
sce_str <- readRDS(file = sce_path_inp[[2]])

#-------------------------------------------------------------------------------
# separate objects

hsc_subclusters <- c(1:length(unique(sce_hsc$annotation_subcluster)))
str_subclusters <- c(1:length(unique(sce_str$annotation_subcluster)))

hsc_subclusters_list <- as.list(hsc_subclusters)
names(hsc_subclusters_list) <- hsc_subclusters
str_subclusters_list <- as.list(str_subclusters)
names(str_subclusters_list) <- str_subclusters

#-------------------------------------------------------------------------------
# save separated objects

path_out_hsc <- path_out[grep("hsc", path_out)]
path_out_str <- path_out[grep("str", path_out)]

names(path_out_hsc) <- c(1:length(path_out_hsc))
print(path_out_hsc)
for(i in 1:length(hsc_subclusters)){
  output_path <- path_out_hsc[i]
  print(output_path)
  saveRDS(hsc_subclusters_list[[i]], file = output_path)
}

names(path_out_str) <- c(1:length(path_out_str))
print(path_out_str)
for(i in 1:length(str_subclusters)){
  output_path <- path_out_str[i]
  print(output_path)
  saveRDS(str_subclusters_list[[i]], file = output_path)
}

sessionInfo()
