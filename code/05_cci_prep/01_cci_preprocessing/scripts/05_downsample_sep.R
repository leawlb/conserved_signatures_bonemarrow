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
  name <- paste0(sce_list[[i]]$Species_ID[1], "_", sce_list[[i]]$Age_ID[1])
  sce_list[[i]]$name <- name
  
  names(sce_list)[i] <- sce_list[[i]]$name[1]
}


for(i in 1:length(sce_list)){
  name <- paste0(sce_list[[i]]$Species_ID[1], "_", sce_list[[i]]$Age_ID[1])
  sce_list[[i]]$name <- name
  
  names(sce_list)[i] <- sce_list[[i]]$name[1]
}

print(sce_list)

sce_list_sep <- lapply(sce_list, function(sce){
  name <- sce$name[1]
  
  sce_emi <- sce[,sce$Assignment == "emitter"]
  sce_rec <- sce[,sce$Assignment == "receiver"]
  
  sce_emi$temp_ident <- paste0(name, "_emi")
  sce_rec$temp_ident <- paste0(name, "_rec")
  
  return(list(sce_emi, sce_rec))
})

sce_list_sep <- unlist(sce_list_sep)

for(i in 1:length(sce_list_sep)){
  names(sce_list_sep)[i] <- sce_list_sep[[i]]$temp_ident[1]
}

down_output_list <- downsampleBatches(sce_list_sep, assay.type = "counts")

for(i in 1:length(sce_list_sep)){
  assays(sce_list_sep[[i]])$downsampled <- down_output_list[[i]]
}
print(sce_list_sep[[i]])

for(i in 1:length(sce_list)){
  
  name <- sce_list[[i]]$name[1]
  print(name)
  print(sce_output_path)
  sce_down <- cbind(sce_list_sep[[paste0(name, "_emi")]], 
                    sce_list_sep[[paste0(name, "_rec")]])
  
  stopifnot(c(grepl(sce_down$Age_ID[1], sce_output_path[[i]]),
              grepl(sce_down$Species_ID[1], sce_output_path[[i]])))
  
  saveRDS(sce_down, file = sce_output_path[[i]])
}  
