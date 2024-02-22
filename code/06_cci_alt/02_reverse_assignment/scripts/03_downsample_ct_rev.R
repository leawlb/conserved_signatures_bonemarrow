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
emitters <- unique(assignment_list$Identity[assignment_list$Assignment == "emitter"])
receivers <- unique(assignment_list$Identity[assignment_list$Assignment == "receiver"])
print(emitters)
print(receivers)

# now separate each SCE into different cell types, based on Assignment
sce_list_ct_sep <- lapply(sce_list, function(sce){

  emi_list <- as.list(emitters)
  rec_list <- as.list(receivers)
  
  sce_ct_list_emi <- lapply(emi_list, function(emi){
    if(emi %in% sce$Identity){
      sce_ct <- sce[,which(sce$Identity == emi)]
    return(sce_ct)
    }
  })
  sce_ct_list_rec <- lapply(rec_list, function(rec){
    print(rec)
    if(rec %in% sce$Identity){
      print("yes")
      sce_ct <- sce[,which(sce$Identity == rec)]
      return(sce_ct)
    }
  })
  names(sce_ct_list_emi) <- unlist(emi_list)
  names(sce_ct_list_rec) <- unlist(rec_list)
  
  # remove list items that are NULL because these cell types did not exist
  sce_ct_list_emi <- Filter(Negate(is.null), sce_ct_list_emi)
  sce_ct_list_rec <- Filter(Negate(is.null), sce_ct_list_rec)
  
  return(list("emitters" = sce_ct_list_emi, "receivers" = sce_ct_list_rec))
})

# sort it into one list containing all emitters, and one containing all receivers
conds <- names(sce_list)
sce_list_ct_sep_sorted <- list()
for(c in conds){
  for(e in emitters){
    if(e %in% names(sce_list_ct_sep[[c]]$emitters)){
    sce_list_ct_sep_sorted[["emitters"]][[e]][[c]] <- sce_list_ct_sep[[c]]$emitters[[e]]
    }
  }
  for(r in receivers){
    if(r %in% names(sce_list_ct_sep[[c]]$receivers)){
    sce_list_ct_sep_sorted[["receivers"]][[r]][[c]] <- sce_list_ct_sep[[c]]$receivers[[r]]
    }
  }  
}

sce_list_ct_sep_emi <- unlist(sce_list_ct_sep_sorted$emitters)
sce_list_ct_sep_rec <- unlist(sce_list_ct_sep_sorted$receivers)
print(names(sce_list_ct_sep_emi))
print(names(sce_list_ct_sep_rec))

# downsample each list
down_output_list_emi <- downsampleBatches(sce_list_ct_sep_emi, assay.type = "counts")
down_output_list_rec <- downsampleBatches(sce_list_ct_sep_rec, assay.type = "counts")
stopifnot(!identical(down_output_list_emi[[1]], counts(sce_list_ct_sep_emi[[1]])))

for(i in 1:length(sce_list_ct_sep_emi)){
  assays(sce_list_ct_sep_emi[[i]])$downsampled <- down_output_list_emi[[i]]
}
for(i in 1:length(sce_list_ct_sep_rec)){
  assays(sce_list_ct_sep_rec[[i]])$downsampled <- down_output_list_rec[[i]]
}

print(sce_list_ct_sep_emi)
print(sce_list_ct_sep_emi[[1]])

# cbind the sce_lists together back to the original SCE but now with downsampled counts
for(i in 1:length(sce_list)){
  
  name <- sce_list[[i]]$name[1]
  print(name)
  #print(sce_output_path)
  
  sce_ct_list_emi <- sce_list_ct_sep_emi[grep(name, names(sce_list_ct_sep_emi))]
  sce_ct_list_rec <- sce_list_ct_sep_rec[grep(name, names(sce_list_ct_sep_rec))]
  print(sce_ct_list_emi)
  # make empty SCE to start cbind
  sce_down <- sce_ct_list_emi[[1]][,0]
  
  for(j in 1:length(sce_ct_list_emi)){
    sce_down <- cbind(sce_down, sce_ct_list_emi[[j]])
  }
  for(j in 1:length(sce_ct_list_rec)){
    sce_down <- cbind(sce_down, sce_ct_list_rec[[j]])
  }
  
  # make sure that the objects are otherwise identical, just with downsampled
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
