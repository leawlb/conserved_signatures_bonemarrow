#-------------------------------------------------------------------------------

library(stringr)
library(DropletUtils)
set.seed(37)
source(file = "../../source/sce_functions.R")

sce_05_path <- snakemake@input[["sce_05_path"]] # cell annotation output
individuals <- snakemake@params[["individuals"]] # all objects to merge
samples_to_remove <- snakemake@params[["samples_to_remove"]] # samples to remove from merging

species_curr <- snakemake@wildcards[["species"]]
print(species_curr)

# only keep objects that are not to be removed
individuals <- individuals[!individuals %in% samples_to_remove]
individuals_curr <- individuals[grepl(species_curr, individuals)]
print(individuals_curr)

# load SCE objects
print(sce_05_path)

sce_list <- list()

for(i in individuals_curr){
  print(paste0(sce_05_path, "/sce_", i, "-05"))
      sce_list[[i]] <- readRDS(file = paste0(sce_05_path, "/sce_", i, "-05"))
}

# only keep objects that are not to be removed (failsave)
sce_list <- sce_list[names(sce_list)[!names(sce_list) %in% samples_to_remove]]

#-------------------------------------------------------------------------------
## prepare fraction-wise merge (within species)

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

sce_list_old <- sce_list[grep("old", names(sce_list))]
sce_list_yng <- sce_list[grep("yng", names(sce_list))]

print(names(sce_list_old))
print(names(sce_list_yng))

# discard designated objects as specified in metadata
# only merge all objects if specified in config

merge_lists <- function(sce_list){

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
  print(ncol(sce_merged))
  return(sce_merged)
}

sce_old_merged <- merge_lists(sce_list_old)
sce_yng_merged <- merge_lists(sce_list_yng)

#-------------------------------------------------------------------------------
# dimensionality reduction

set.seed(37)
sce_old_merged <- reduce_dims(sce_old_merged,
                              nr_hvgs = snakemake@params[["nr_hvgs"]])
sce_yng_merged <- reduce_dims(sce_yng_merged,
                              nr_hvgs = snakemake@params[["nr_hvgs"]])

test_merged <- function(sce_merged){
  print(unique(sce_merged$Age_ID))
  print(unique(sce_merged$Species_ID))
  print(unique(sce_merged$Object_ID))
}

test_merged(sce_old_merged)
test_merged(sce_yng_merged)

saveRDS(sce_old_merged, file = snakemake@output[["sce_06_old"]]) # merging output
saveRDS(sce_yng_merged, file = snakemake@output[["sce_06_yng"]]) # merging output
