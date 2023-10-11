library(SingleCellExperiment)
library(scuttle)

set.seed(37)

#-------------------------------------------------------------------------------

sce_input_path <- snakemake@input[["sce_input_path"]]
sce_output_path <- snakemake@output[["sce_output_path"]]
#print(sce_input_path)
#print(sce_output_path)

sce_list <- list()

for(i in 1:length(sce_input_path)){
  print(sce_input_path[[i]])
  sce_list[[i]] <- readRDS(file = sce_input_path[[i]])
  name <- paste0(sce_list[[i]]$Species_ID[1], "_", sce_list[[i]]$Age_ID[1])
  print(name)
  sce_list[[i]]$name <- name
  
  names(sce_list)[i] <- sce_list[[i]]$name[1]
}

print(sce_list)


for(i in 1:length(sce_list)){
  name <- paste0(sce_list[[i]]$Species_ID[1], "_", sce_list[[i]]$Age_ID[1])
  sce_list[[i]]$name <- name
  names(sce_list)[i] <- name
}

print(sce_list)

# now separate into different cell types
sce_list_ct_sep <- lapply(sce_list, function(sce){
  print(sce$name[1])
  cts_list <- as.list(unique(sce$Identity))
  
  sce_ct_list <- lapply(cts_list, function(cts){
    sce_ct <- sce[,which(sce$Identity == cts)]
    return(sce_ct)
  })
  names(sce_ct_list) <- unlist(cts_list)
  return(sce_ct_list)
})


sce_list_ct_sep <- unlist(sce_list_ct_sep)
print(sce_list_ct_sep)

down_output_list <- downsampleBatches(sce_list_ct_sep, assay.type = "counts")

for(i in 1:length(sce_list_ct_sep)){
  assays(sce_list_ct_sep[[i]])$downsampled <- down_output_list[[i]]
}
print(sce_list_ct_sep[[i]])

for(i in 1:length(sce_list)){
  
  name <- sce_list[[i]]$name[1]
  print(name)
  print(sce_output_path)
  sce_ct_list <- sce_list_ct_sep[grep(name, names(sce_list_ct_sep))]
  sce_down <- sce_ct_list[[1]][,0]
  
  for(j in 1:length(sce_ct_list)){
    sce_down <- cbind(sce_down, sce_ct_list[[j]])
  }
  print(dim(sce_down))
  print(dim(sce_list[[i]]))
  
  stopifnot(c(grepl(sce_down$Age_ID[1], sce_output_path[[i]]),
              grepl(sce_down$Species_ID[1], sce_output_path[[i]])))
  
  saveRDS(sce_down, file = sce_output_path[[i]])
}  
