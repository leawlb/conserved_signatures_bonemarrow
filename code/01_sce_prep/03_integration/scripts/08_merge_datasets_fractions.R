#-------------------------------------------------------------------------------

library(stringr)
library(scran)
set.seed(37)

source(file = snakemake@params[["functions"]])

#-------------------------------------------------------------------------------

sce_input_path <- snakemake@input[["sce_input_path"]] 

# all objects IDs/individuals to merge
individuals <- snakemake@params[["individuals"]]
# samples to exclude from merging
samples_to_remove <- snakemake@params[["samples_to_remove"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]

# only keep objects that are not to be removed
individuals <- individuals[!individuals %in% samples_to_remove]
print(individuals)

# iteratively load required SCE objects only
print(sce_input_path)

sce_list <- list()

for(j in 1:length(sce_input_path)){
  # since species is not in a wildcard, it needs to be extracted from the path
  # need to do it per species as in 01/01 
  species_curr <- stringr::str_sub(sce_input_path[j], start= -4) 
  for(i in individuals){
    if(grepl(species_curr, i) == TRUE){
      sce_list[[i]] <- base::readRDS(file = base::paste0(sce_input_path[j], 
                                                         "/sce_", 
                                                         i,
                                                         "-07"))
    }
  }
}
# only keep objects that are not to be removed (failsave)
sce_list <- sce_list[names(sce_list)[!names(sce_list) %in% samples_to_remove]]

#-------------------------------------------------------------------------------
## prepare fraction-wise merge

# get genes shared by all objects
rownames_list <- base::invisible(lapply(sce_list, function(sce){
  print(sce$Object_ID[1]) # stop lapply from printing rownames
  rownames <- rowData(sce)$Symbol
  return(rownames)
}))

print(rownames_list[[1]][1:10])

shared_genes <- BiocGenerics::Reduce(intersect, rownames_list)
print(length(shared_genes))

# subset all objects by shared_genes
sce_list <- lapply(sce_list, function(sce){
  sce <- sce[rownames(sce) %in% shared_genes,]
  return(sce)
})

# add individual sample name to barcodes to avoid random doublets
sce_list <- lapply(sce_list, function(sce){
  colnames(sce) <- base::paste0(colnames(sce), "_",  sce$Object_ID)
  return(sce)
})

#-------------------------------------------------------------------------------
## fraction merge  

sce_list_hsc <- sce_list[base::grep("hsc", names(sce_list))]
sce_list_str <- sce_list[base::grep("str", names(sce_list))]

print("HSC")
print(sce_list_hsc)
print("STR")
print(sce_list_str)

# only merge objects if specified in config
for(i in 1:length(sce_list_hsc)){
  if(i == 1){
    if(sce_list_hsc[[i]]$Keep_sample[1] == TRUE){ 
        sce_merged_hsc <- sce_list_hsc[[i]] # set first item of list
      }else{ # only if not to be discarded
        stop("first item in list is to be discarded")
      }
    }else{
      if(sce_list_hsc[[i]]$Keep_sample[1] == TRUE){
        sce_merged_hsc <- BiocGenerics::cbind(sce_merged_hsc, sce_list_hsc[[i]]) 
    }
  }
}
print(nrow(sce_merged_hsc))

for(i in 1:length(sce_list_str)){
  if(i == 1){
    if(sce_list_str[[i]]$Keep_sample[1] == TRUE){ 
      sce_merged_str <- sce_list_str[[i]] 
    }else{ 
      stop("first item in list is to be discarded")
    }
  }else{
    if(sce_list_str[[i]]$Keep_sample[1] == TRUE){
      sce_merged_str <- BiocGenerics::cbind(sce_merged_str, sce_list_str[[i]]) 
    }
  }
}
print(nrow(sce_merged_str))

#-------------------------------------------------------------------------------
# dimensionality reduction (own function)

sce_merged_hsc <- reduce_dims(sce_merged_hsc, nr_hvgs = nr_hvgs)
sce_merged_str <- reduce_dims(sce_merged_str, nr_hvgs = nr_hvgs)

base::saveRDS(sce_merged_hsc, file = snakemake@output[["sce_output_hsc"]]) 
base::saveRDS(sce_merged_str, file = snakemake@output[["sce_output_str"]]) 

utils::sessionInfo()