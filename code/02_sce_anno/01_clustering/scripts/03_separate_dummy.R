#-------------------------------------------------------------------------------
# save dummy object to allow using "cluster" as wildcard later on

set.seed(37)

sce_path_inp <- snakemake@input[["sce_input"]]
sce_path_out <- snakemake@output[["output"]]

#-------------------------------------------------------------------------------

sce_hsc <- base::readRDS(file = sce_path_inp[[1]])
sce_str <- base::readRDS(file = sce_path_inp[[2]])

#-------------------------------------------------------------------------------
# separate objects

hsc_clusters <- c(1:length(base::unique(sce_hsc$cluster_louvain)))
str_clusters <- c(1:length(base::unique(sce_str$cluster_louvain)))

hsc_clusters_list <- as.list(hsc_clusters)
str_clusters_list <- as.list(str_clusters)

names(hsc_clusters_list) <- hsc_clusters
names(str_clusters_list) <- str_clusters

#-------------------------------------------------------------------------------
# save in place of SCE objects

sce_path_out_hsc <- sce_path_out[base::grep("hsc", sce_path_out)]
sce_path_out_str <- sce_path_out[base::grep("str", sce_path_out)]

names(sce_path_out_hsc) <- c(1:length(sce_path_out_hsc))

print("hsc")
print(sce_path_out_hsc)
for(i in hsc_clusters){
  print(i)
  output_path <- sce_path_out_hsc[i]
  print(output_path)
  base::saveRDS(hsc_clusters_list[[i]], file = output_path)
}

print("str")
print(sce_path_out_str)
for(i in str_clusters){
  print(i)
  output_path <- sce_path_out_str[i]
  print(output_path)
  base::saveRDS(str_clusters_list[[i]], file = output_path)
}

utils::sessionInfo()
