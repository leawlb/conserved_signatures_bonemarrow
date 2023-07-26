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
# separate objects

sce_hsc$annotation_subcluster <- factor(sce_hsc$annotation_subcluster, 
                                        levels = unique(sce_hsc$annotation_subcluster))
sce_str$annotation_subcluster <- factor(sce_str$annotation_subcluster, 
                                        levels = unique(sce_str$annotation_subcluster))

print(levels(sce_hsc$annotation_subcluster))
print(levels(sce_str$annotation_subcluster))

sce_hsc_list <- list()
for(i in 1:length(levels(sce_hsc$annotation_subcluster))){
  sce_hsc_list[[i]] <- sce_hsc[,sce_hsc$annotation_subcluster == levels(
    sce_hsc$annotation_subcluster)[i]]
}

sce_str_list <- list()
for(i in 1:length(levels(sce_str$annotation_subcluster))){
  sce_str_list[[i]] <- sce_str[,sce_str$annotation_subcluster == levels(
    sce_str$annotation_subcluster)[i]]
}

#-------------------------------------------------------------------------------
# save separated objects

sce_path_out_hsc <- sce_path_out[grep("hsc", sce_path_out)]
sce_path_out_str <- sce_path_out[grep("str", sce_path_out)]

names(sce_path_out_hsc) <- c(1:length(sce_path_out_hsc))

print(sce_path_out_hsc)
print(levels(sce_hsc$annotation_subcluster))
for(i in 1:length(levels(sce_hsc$annotation_subcluster))){
  print(levels(sce_hsc$annotation_subcluster)[i])
  output_path <- sce_path_out_hsc[i]
  print(output_path)
  saveRDS(sce_hsc_list[[i]], file = output_path)
}
print("finished HSCs")

names(sce_path_out_str) <- c(1:length(sce_path_out_str))
for(i in 1:length(levels(sce_str$annotation_subcluster))){
  print(levels(sce_str$annotation_subcluster)[i])
  output_path <- sce_path_out_str[i]
  print(output_path)
  saveRDS(sce_str_list[[i]], file = output_path)
}
