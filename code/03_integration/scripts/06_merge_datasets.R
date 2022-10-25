#-------------------------------------------------------------------------------

library(stringr)

sce_04_path <- snakemake@input[["sce_04"]] # correct preprocessing output
sce_06_path <- snakemake@output[["sce_06"]] # correct integration output
individuals <- snakemake@params[["individuals"]] # all objects to merge

species_curr <- snakemake@wildcards[["species"]]

sce_list_04 <- list()

if(!is.null(species_curr)){ # if species wildcard was used (merge_species)
  for(i in individuals){
    if(grepl(species_curr, i) == TRUE){
      sce_list_04[[i]] <- readRDS(file = paste0(sce_04_path, "/sce_", i, "-04"))
    }
  }
}

if(is.null(species_curr)){ # if no wildcard was used (merge all)
  for(j in 1:length(sce_04_path)){
    # the "current" species i is the last part of the path
    species_curr <- str_sub(sce_04_path[j], start= -4) 
    for(i in individuals){
      if(grepl(species_curr, i) == TRUE){
        sce_list_04[[i]] <- readRDS(file = paste0(sce_04_path[j], 
                                                  "/sce_", i, "-04"))
      }
    }
  }
}

print(names(sce_list_04))

#-------------------------------------------------------------------------------
## prepare merge

# get genes shared by all objects
rownames_list <- lapply(sce_list_04, function(sce){
  
  print(rownames(sce)[1]) # this only works with this line of code here? why? DO NOT DELETE 
  rn <- rownames(sce)
  return(rn)
})
print(rownames_list)[[1]][1:10] 
print(names(rownames_list))
shared_genes <- Reduce(intersect, rownames_list)

# subset all objects by shared_genes
sce_list_04 <- lapply(sce_list_04, function(sce){
  sce <- sce[rownames(sce) %in% shared_genes,]
  return(sce)
})

# add individual sample name to barcodes to avoid random doublets
sce_list_04 <- lapply(sce_list_04, function(sce){
  colnames(sce) <- paste0(colnames(sce), "_", sce$Object_ID)
  return(sce)
})

#-------------------------------------------------------------------------------
## all merge  

# discard designated objects as specified in metadata
# only merge all objects if specified in config
for(i in 1:length(sce_list_04)){
  if(i == 1){
    if(sce_list_04[[i]]$Keep_sample[1] == TRUE){ 
        sce_merged <- sce_list_04[[i]] # set first item of list
      }else{ # only if not to be discarded
        stop("first item in list is to be discarded")
      }
    }else{
      if(sce_list_04[[i]]$Keep_sample[1] == TRUE){
        sce_merged <- cbind(sce_merged, sce_list_04[[i]]) 
    }
  }
}
print(nrow(sce_merged))

saveRDS(sce_merged, file = sce_06_path)