#-------------------------------------------------------------------------------

library(scran)
set.seed(37)

source(file = snakemake@params[["functions"]])

#-------------------------------------------------------------------------------

sce_input_path <- snakemake@input[["sce_input_path"]] 

individuals <- snakemake@params[["individuals"]]
samples_to_remove <- snakemake@params[["samples_to_remove"]]
nr_hvgs <- snakemake@params[["nr_hvgs"]]

species_curr <- snakemake@wildcards[["species"]]
print(species_curr)

individuals <- individuals[!individuals %in% samples_to_remove]
individuals_curr <- individuals[base::grepl(species_curr, individuals)]
print(individuals_curr)

print(sce_input_path)

sce_list <- list()
for(i in individuals_curr){
  print(base::paste0(sce_input_path, 
                     "/sce_",
                     i, 
                     "-07"))
      sce_list[[i]] <- base::readRDS(file = base::paste0(sce_input_path, 
                                                         "/sce_",
                                                         i,
                                                         "-07"))
}

sce_list <- sce_list[names(sce_list)[!names(sce_list) %in% samples_to_remove]]

#-------------------------------------------------------------------------------
## prepare fraction-wise merge (within species)

rownames_list <- base::invisible(lapply(sce_list, function(sce){
  print(sce$Object_ID[1]) # stop lapply from printing rownames
  rownames <- rowData(sce)$Symbol
  return(rownames)
}))

print(rownames_list[[1]][1:10])

shared_genes <- BiocGenerics::Reduce(intersect, rownames_list)
print(length(shared_genes))

sce_list <- lapply(sce_list, function(sce){
  sce <- sce[rownames(sce) %in% shared_genes,]
  return(sce)
})

sce_list <- lapply(sce_list, function(sce){
  colnames(sce) <- base::paste0(colnames(sce), "_", sce$Object_ID)
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

for(i in 1:length(sce_list_hsc)){
  if(i == 1){
    if(sce_list_hsc[[i]]$Keep_sample[1] == TRUE){ 
      sce_merged_hsc <- sce_list_hsc[[i]] 
    }else{ 
      stop("first item in list is to be discarded")
    }
  }else{
    if(sce_list_hsc[[i]]$Keep_sample[1] == TRUE){
      sce_merged_hsc <- BiocGenerics::cbind(sce_merged_hsc, sce_list_hsc[[i]]) 
    }
  }
}
print(nrow(sce_merged_hsc))
print(ncol(sce_merged_hsc))

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
print(ncol(sce_merged_str))

#-------------------------------------------------------------------------------

set.seed(37)
sce_merged_hsc <- reduce_dims(sce_merged_hsc, nr_hvgs = nr_hvgs)
sce_merged_str <- reduce_dims(sce_merged_str, nr_hvgs = nr_hvgs)

base::saveRDS(sce_merged_hsc, file = snakemake@output[["sce_output_hsc"]]) 
base::saveRDS(sce_merged_str, file = snakemake@output[["sce_output_str"]]) 

utils::sessionInfo()