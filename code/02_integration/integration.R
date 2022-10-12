#-------------------------------------------------------------------------------

setwd("~/Interspecies_BM_phd/code/02_integration")
library(config)

sce_04_path <- snakemake@input[["sce_04_path"]] # correct preprocessing output
sce_06_path <- snakemake@output[["sce_06_path"]] # correct integration output
targets <- snakemake@params[["targets"]] # all objects to merge
run_sce_all<- snakemake@params[["run_sce_all"]] # all objects to merge

sce_list_04 <- list()

for(i in targets){
  sce_list_04[[i]] <- readRDS(file = paste0(sce_04_path, "/sce_", i, "-04"))
}

#-------------------------------------------------------------------------------
## prepare merge

# get genes shared by all objects
rownames_list <- lapply(sce_list_04, function(sce){
  rownames <- rownames(sce)
  return(rownames)
})
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
if(run_sce_all == TRUE){
  for(i in 1:length(sce_list_04)){
    if(i == 1){
      if(sce_list_04[[i]]$Discard[1] == TRUE){ 
        stop("first item in list is to be discarded")
      }else{ # only if not to be discarded
        sce_all <- sce_list_04[[i]] # set first item of list
      }
    }else{
      if(sce_list_04[[i]]$Discard[1] == FALSE){
        sce_all <- cbind(sce_all, sce_list_04[[i]]) # successively add items
      }
    }
  }
}

#-------------------------------------------------------------------------------
## species-wise merge 

# it is necessary to skip over species that are not found in this run

# mcar
sce_list_mcar <- sce_list_04[targets[grep("mcar", targets)]]

if(length(sce_list_mcar) > 0){ 
  for(i in 1:length(sce_list_mcar)){
    if(i == 1){ 
      if(sce_list_mcar[[i]]$Discard[1] == TRUE){
        stop("first item in list is to be discarded")
      }else{
        sce_mcar <- sce_list_mcar[[i]]
      }
    }else if(i > 0){
      if(sce_list_mcar[[i]]$Discard[1] == FALSE){
        sce_mcar <- cbind(sce_mcar, sce_list_mcar[[i]])
      }
    }
  }
}

# mcas
sce_list_mcas <- sce_list_04[targets[grep("mcas", targets)]]

if(length(sce_list_mcas) > 0){ 
  for(i in 1:length(sce_list_mcas)){
    if(i == 1){ 
      if(sce_list_mcas[[i]]$Discard[1] == TRUE){
        stop("first item in list is to be discarded")
      }else{
        sce_mcas <- sce_list_mcas[[i]]
      }
    }else if(i > 0){
      if(sce_list_mcas[[i]]$Discard[1] == FALSE){
        sce_mcas <- cbind(sce_mcas, sce_list_mcas[[i]])
      }
    }
  }
}

# mmus
sce_list_mmus <- sce_list_04[targets[grep("mmus", targets)]]

if(length(sce_list_mmus) > 0){ 
  for(i in 1:length(sce_list_mmus)){
    if(i == 1){ 
      if(sce_list_mmus[[i]]$Discard[1] == TRUE){
        stop("first item in list is to be discarded")
      }else{
        sce_mmus <- sce_list_mmus[[i]]
      }
    }else if(i > 0){
      if(sce_list_mmus[[i]]$Discard[1] == FALSE){
        sce_mmus <- cbind(sce_mmus, sce_list_mmus[[i]])
      }
    }
  }
}

# mspr
sce_list_mspr <- sce_list_04[targets[grep("mspr", targets)]]
if(length(sce_list_mspr) > 0){ 
  for(i in 1:length(sce_list_mspr)){
    if(i == 1){ 
      if(sce_list_mspr[[i]]$Discard[1] == TRUE){
        stop("first item in list is to be discarded")
      }else{
        sce_mspr <- sce_list_mspr[[i]]
      }
    }else if(i > 0){
      if(sce_list_mspr[[i]]$Discard[1] == FALSE){
        sce_mspr <- cbind(sce_mspr, sce_list_mspr[[i]])
      }
    }
  }
}

if(exists("sce_mspr")){
  print("yes")
}

#-------------------------------------------------------------------------------
## save under new names which can from now on be used as wildcards
saveRDS(sce_all, file = paste0(sce_06_path, "/sce_all-06"))
if(exists("sce_mcar")){
  saveRDS(sce_mcar, file = paste0(sce_06_path, "/sce_mcar-06"))
}
if(exists("sce_mcas")){
  saveRDS(sce_mcas, file = paste0(sce_06_path, "/sce_mcas-06"))
}
if(exists("sce_mmus")){
  saveRDS(sce_mmus, file = paste0(sce_06_path, "/sce_mmus-06"))
}
if(exists("sce_mspr")){
  saveRDS(sce_mspr, file = paste0(sce_06_path, "/sce_mspr-06"))
}