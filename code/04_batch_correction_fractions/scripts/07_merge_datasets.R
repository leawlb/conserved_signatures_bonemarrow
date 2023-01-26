#-------------------------------------------------------------------------------

library(stringr)
library(DropletUtils)
set.seed(37)
source(file = snakemake@params[["sce_functions"]])

sce_06_path <- snakemake@input[["sce_06"]] # correct cell annotation output
sce_07_hsc_path <- snakemake@output[["sce_07_hsc"]] # correct merging output
sce_07_str_path <- snakemake@output[["sce_07_str"]] # correct merging output
individuals <- snakemake@params[["individuals"]] # all objects to merge
samples_to_remove <- snakemake@params[["samples_to_remove"]] # samples to remove from merging
#print(samples_to_remove)  
#print(individuals)

# only keep objects that are not to be removed
individuals <- individuals[!individuals %in% samples_to_remove]
print(individuals)

# load SCE objects
print(sce_06_path)

sce_list <- list()

for(j in 1:length(sce_06_path)){
  species_curr <- str_sub(sce_06_path[j], start= -4) 
  for(i in individuals){
    if(grepl(species_curr, i) == TRUE){
      sce_list[[i]] <- readRDS(file = paste0(sce_06_path[j], 
                                             "/sce_", i, "-06"))
    }
  }
}

# only keep objects that are not to be removed (failsave)
sce_list <- sce_list[names(sce_list)[!names(sce_list) %in% samples_to_remove]]

#-------------------------------------------------------------------------------
## prepare fraction-wise merge

# get genes shared by all objects
rownames_list <- invisible(lapply(sce_list, function(sce){
  print(sce$Object_ID[1]) # stop lapply from printing rownames
  rownames <- rowData(sce)$Symbol
  return(rownames)
}))

print(rownames_list[[1]][1:10])

shared_genes <- Reduce(intersect, rownames_list)

# subset all objects by shared_genes
sce_list <- lapply(sce_list, function(sce){
  sce <- sce[rownames(sce) %in% shared_genes,]
  return(sce)
})

# add individual sample name to barcodes to avoid random doublets
sce_list <- lapply(sce_list, function(sce){
  colnames(sce) <- paste0(colnames(sce), "_", sce$Object_ID)
  return(sce)
})

#-------------------------------------------------------------------------------
## fraction merge  

sce_list_hsc <- sce_list[grep("hsc", names(sce_list))]
sce_list_str <- sce_list[grep("str", names(sce_list))]

print("HSC")
print(sce_list_hsc)
print("STR")
print(sce_list_str)

# discard designated objects as specified in metadata
# only merge all objects if specified in config
for(i in 1:length(sce_list_hsc)){
  if(i == 1){
    if(sce_list_hsc[[i]]$Keep_sample[1] == TRUE){ 
        sce_merged_hsc <- sce_list_hsc[[i]] # set first item of list
      }else{ # only if not to be discarded
        stop("first item in list is to be discarded")
      }
    }else{
      if(sce_list_hsc[[i]]$Keep_sample[1] == TRUE){
        sce_merged_hsc <- cbind(sce_merged_hsc, sce_list_hsc[[i]]) 
    }
  }
}
print(nrow(sce_merged_hsc))

for(i in 1:length(sce_list_str)){
  if(i == 1){
    if(sce_list_str[[i]]$Keep_sample[1] == TRUE){ 
      sce_merged_str <- sce_list_str[[i]] # set first item of list
    }else{ # only if not to be discarded
      stop("first item in list is to be discarded")
    }
  }else{
    if(sce_list_str[[i]]$Keep_sample[1] == TRUE){
      sce_merged_str <- cbind(sce_merged_str, sce_list_str[[i]]) 
    }
  }
}
print(nrow(sce_merged_str))

#-------------------------------------------------------------------------------
# dimensionality reduction

set.seed(37)
sce_merged_hsc <- reduce_dims(sce_merged_hsc,
                              nr_hvgs = snakemake@params[["nr_hvgs"]])
sce_merged_str <- reduce_dims(sce_merged_str,
                              nr_hvgs = snakemake@params[["nr_hvgs"]])

saveRDS(sce_merged_hsc, file = sce_07_hsc_path)
saveRDS(sce_merged_str, file = sce_07_str_path)