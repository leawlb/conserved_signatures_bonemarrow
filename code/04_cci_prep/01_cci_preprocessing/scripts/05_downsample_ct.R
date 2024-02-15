library(SingleCellExperiment)
library(scuttle)

set.seed(37)

#-------------------------------------------------------------------------------

sce_input_path <- snakemake@input[["sce_input_path"]]
sce_output_path <- snakemake@output[["sce_output_path"]]
assignment_list <- read.csv(snakemake@input[["assignment"]], 
                            header = TRUE, 
                            sep = ";", 
                            check.names=FALSE, 
                            stringsAsFactors=FALSE, 
                            as.is=TRUE, 
                            colClasses = "character")

sce_list <- list()

# load list of SCE objects to downsample
for(i in 1:length(sce_input_path)){
  print(sce_input_path[[i]])
  sce_list[[i]] <- readRDS(file = sce_input_path[[i]])
  name <- paste0(sce_list[[i]]$Species_ID[1], "_", sce_list[[i]]$Age_ID[1])
  print(name)
  sce_list[[i]]$name <- name
  
  names(sce_list)[i] <- sce_list[[i]]$name[1]
}

# add name to each list item
for(i in 1:length(sce_list)){
  name <- paste0(sce_list[[i]]$Species_ID[1], "_", sce_list[[i]]$Age_ID[1])
  sce_list[[i]]$name <- name
}
print(names(sce_list))

# extract emitter and receiver cell types to downsample them separately
idents <- unique(assignment_list$Identity)
print(idents)

# now separate each SCE into different cell types, based on Assignment
sce_list_ct_sep <- lapply(sce_list, function(sce){

  ident_list <- as.list(idents)

  sce_ct_list <- lapply(ident_list, function(ct){
    if(ct %in% sce$Identity){
      sce_ct <- sce[,which(sce$Identity == ct)]
    return(sce_ct)
    }
  })
  
  names(sce_ct_list) <- unlist(ident_list)

  # remove list items that are NULL because these cell types did not exist
  sce_ct_list <- Filter(Negate(is.null), sce_ct_list)

  return(sce_ct_list)
})

sce_list_ct_sep <- unlist(sce_list_ct_sep)
print(names(sce_list_ct_sep))
print(sce_list_ct_sep[[1]])

# downsample each list
down_output_list <- downsampleBatches(sce_list_ct_sep, assay.type = "counts")
stopifnot(!identical(down_output_list[[1]], counts(sce_list_ct_sep[[1]])))

for(i in 1:length(sce_list_ct_sep)){
  assays(sce_list_ct_sep[[i]])$downsampled <- down_output_list[[i]]
}

print(names(sce_list_ct_sep))
print(sce_list_ct_sep[[1]])

# cbind the sce_lists together back to the original SCE but now with downsampled counts
for(i in 1:length(sce_list)){
  
  name <- sce_list[[i]]$name[1]
  print(name)
  #print(sce_output_path)
  
  sce_ct_list <- sce_list_ct_sep[grep(name, names(sce_list_ct_sep))]
  
  print(sce_ct_list)
  # make empty SCE to start cbind
  sce_down <- sce_ct_list[[1]][,0]
  
  for(j in 1:length(sce_ct_list)){
    sce_down <- cbind(sce_down, sce_ct_list[[j]])
  }
  
  # make sure that the objects are otherwise identical, just with "downsampled"
  # assay added
  sce_down <- sce_down[,match(colnames(sce_list[[i]]), colnames(sce_down))]
  print(assays(sce_down))

  print(dim(sce_down))
  print(dim(sce_list[[i]]))
  stopifnot(identical(colData(sce_list[[i]]), colData(sce_down)))
  stopifnot(identical(counts(sce_list[[i]]), counts(sce_down)))

  stopifnot(c(grepl(sce_down$Age_ID[1], sce_output_path[[i]]),
              grepl(sce_down$Species_ID[1], sce_output_path[[i]])))
  
  saveRDS(sce_down, file = sce_output_path[[i]])
}  
