#-------------------------------------------------------------------------------

library(scran)
RNGkind("L'Ecuyer-CMRG") 
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
  print(paste0(sce_input_path, "/sce_", i, "-07"))
      sce_list[[i]] <- base::readRDS(file = base::paste0(sce_input_path,
                                                         "/sce_", 
                                                         i,
                                                         "-07"))
}

sce_list <- sce_list[names(sce_list)[!names(sce_list) %in% samples_to_remove]]

#-------------------------------------------------------------------------------
## prepare fraction-wise merge (within ages)

rownames_list <- base::invisible(lapply(sce_list, function(sce){
  print(sce$Object_ID[1]) 
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
## age merge  

sce_list_old <- sce_list[base::grep("old", names(sce_list))]
sce_list_yng <- sce_list[base::grep("yng", names(sce_list))]

print(names(sce_list_old))
print(names(sce_list_yng))

merge_lists <- function(sce_list){

  for(i in 1:length(sce_list)){
    if(i == 1){
      if(sce_list[[i]]$Keep_sample[1] == TRUE){ 
        sce_merged <- sce_list[[i]]
      }else{ 
        stop("first item in list is to be discarded")
      }
    }else{
      if(sce_list[[i]]$Keep_sample[1] == TRUE){
        sce_merged <- BiocGenerics::cbind(sce_merged, sce_list[[i]]) 
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

set.seed(37)
sce_old_merged <- reduce_dims(sce_old_merged, nr_hvgs = nr_hvgs)
sce_yng_merged <- reduce_dims(sce_yng_merged, nr_hvgs = nr_hvgs)

test_merged <- function(sce_merged){
  print(unique(sce_merged$Age_ID))
  print(unique(sce_merged$Species_ID))
  print(unique(sce_merged$Object_ID))
}

test_merged(sce_old_merged)
test_merged(sce_yng_merged)

base::saveRDS(sce_old_merged, file = snakemake@output[["sce_output_old"]]) 
base::saveRDS(sce_yng_merged, file = snakemake@output[["sce_output_yng"]]) 

utils::sessionInfo()