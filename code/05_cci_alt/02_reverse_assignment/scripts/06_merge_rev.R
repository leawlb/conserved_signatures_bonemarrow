
set.seed(37)

#-------------------------------------------------------------------------------

cci_input_path <- snakemake@input[["cci_input_path"]] # cci prep output

# load CCI objects
print(cci_input_path)

cci_list <- list()

for(i in 1:length(cci_input_path)){
  cci_list[[i]] <- readRDS(file = cci_input_path[[i]])
}
print(cci_list)

#-------------------------------------------------------------------------------
## prepare fraction-wise merge (within species)

cci_list <- lapply(cci_list, function(cci){
  
  name <- paste0(cci$Identities$Species[1], "_", cci$Identities$Age[1])
  
  colnames(cci$Score) <- paste0(colnames(cci$Score), "_", name)
  colnames(cci$Receptorrank) <- paste0(colnames(cci$Receptorrank), "_", name)
  colnames(cci$Ligandrank) <- paste0(colnames(cci$Ligandrank), "_", name)
  
  
  cci$Identities$name <- rep(name, nrow(cci$Identities))
  return(cci)
})


cci_list_new <- cci_list[[1]]
for(i in 2:length(cci_list)){
  cci_list_new$Score <- cbind(cci_list_new$Score, cci_list[[i]]$Score)
  cci_list_new$Receptorrank <- cbind(cci_list_new$Receptorrank, 
                                      cci_list[[i]]$Receptorrank)
  cci_list_new$Ligandrank <- cbind(cci_list_new$Ligandrank, 
                                    cci_list[[i]]$Ligandrank)
  
  cci_list_new$Identities <- rbind(cci_list_new$Identities, 
                                    cci_list[[i]]$Identities)
  
}

print(cci_list_new$Identities)

saveRDS(cci_list_new, file = snakemake@output[["cci_output"]])
