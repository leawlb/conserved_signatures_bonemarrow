#-------------------------------------------------------------------------------

library(stringr)
library(DropletUtils)

sce_06_path <- snakemake@input[["sce_06"]] # correct cell annotation output
sce_07_path <- snakemake@output[["sce_07"]] # correct merging output
individuals <- snakemake@params[["individuals"]] # all objects to merge
samples_to_remove <- snakemake@params[["samples_to_remove"]] # samples to remove from merging
print(samples_to_remove)  

species_curr <- snakemake@wildcards[["species"]]
source(file = snakemake@params[["sce_functions"]])

sce_list <- list()

if(!is.null(species_curr)){ # if species wildcard was used (merge_species)
  print(species_curr)
  for(i in individuals){
    if(grepl(species_curr, i) == TRUE){
      sce_list[[i]] <- readRDS(file = paste0(sce_06_path, "/sce_", i, "-06"))
    }
  }
}

if(is.null(species_curr)){ # if no wildcard was used (merge all)
  print("all")
  for(j in 1:length(sce_06_path)){
    # the "current" species i is the last part of the path
    species_curr <- str_sub(sce_06_path[j], start= -4) 
    for(i in individuals){
      if(grepl(species_curr, i) == TRUE){
        sce_list[[i]] <- readRDS(file = paste0(sce_06_path[j], 
                                                  "/sce_", i, "-06"))
      }
    }
  }
}

# only keep objects that are not to be removed
#print(names(sce_list)[names(sce_list) %in% samples_to_remove])
sce_list <- sce_list[names(sce_list)[!names(sce_list) %in% samples_to_remove]]

#-------------------------------------------------------------------------------
## prepare merge

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
## all merge  

# discard designated objects as specified in metadata
# only merge all objects if specified in config
for(i in 1:length(sce_list)){
  if(i == 1){
    if(sce_list[[i]]$Keep_sample[1] == TRUE){ 
        sce_merged <- sce_list[[i]] # set first item of list
      }else{ # only if not to be discarded
        stop("first item in list is to be discarded")
      }
    }else{
      if(sce_list[[i]]$Keep_sample[1] == TRUE){
        sce_merged <- cbind(sce_merged, sce_list[[i]]) 
    }
  }
}
print(nrow(sce_merged))

#-------------------------------------------------------------------------------
# dimensionality reduction

sce_merged <- reduce_dims(sce_merged, nr_hvgs = snakemake@params[["nr_hvgs"]])

saveRDS(sce_merged, file = sce_07_path)