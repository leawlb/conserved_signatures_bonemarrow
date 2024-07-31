#-------------------------------------------------------------------------------
# this is just another dummy rule/script to be able to use the new cell types
# as wildcards down the stream

set.seed(37)

#-------------------------------------------------------------------------------
# read objects

sce_path_inp <- snakemake@input[["sce_input"]]
path_out <- snakemake@output[["output"]]

sce_hsc <- base::readRDS(file = sce_path_inp[[1]])
sce_str <- base::readRDS(file = sce_path_inp[[2]])

#-------------------------------------------------------------------------------
# separate objects

# cluster_after_sub contains numeric clusters = number of cell types
# instead of the original louvain clusters, some of which were removed
hsc_subclusters <- c(1:length(base::unique(sce_hsc$cluster_after_sub)))
str_subclusters <- c(1:length(base::unique(sce_str$cluster_after_sub)))

hsc_subclusters_list <- as.list(hsc_subclusters)
str_subclusters_list <- as.list(str_subclusters)

names(hsc_subclusters_list) <- hsc_subclusters
names(str_subclusters_list) <- str_subclusters

print(names(hsc_subclusters_list))
print(names(str_subclusters_list))

#-------------------------------------------------------------------------------
# save separated objects (which is just the cluster nr (after subclustering))

path_out_hsc <- path_out[base::grep("hsc", path_out)]
path_out_str <- path_out[base::grep("str", path_out)]

names(path_out_hsc) <- c(1:length(path_out_hsc))
print(path_out_hsc)
for(i in 1:length(hsc_subclusters)){
  output_path <- path_out_hsc[i]
  print(output_path)
  base::saveRDS(hsc_subclusters_list[[i]], file = output_path)
}

names(path_out_str) <- c(1:length(path_out_str))
print(path_out_str)
for(i in 1:length(str_subclusters)){
  output_path <- path_out_str[i]
  print(output_path)
  base::saveRDS(str_subclusters_list[[i]], file = output_path)
}

utils::sessionInfo()
