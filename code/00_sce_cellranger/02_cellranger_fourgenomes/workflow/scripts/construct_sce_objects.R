#-------------------------------------------------------------------------------

library(tidyverse, quietly = TRUE) 
library(DropletUtils, quietly = TRUE) 
set.seed(1)

individual_curr <- snakemake@wildcards[["individual"]] # currently loaded individual sample 
IDENTIFIERS <- snakemake@params[["identifiers"]] 

# construct SCE object from raw cellranger output
print(paste0(snakemake@input[["output_cellranger"]], "/raw_feature_bc_matrix"))
sce <- read10xCounts(samples = paste0(snakemake@input[["output_cellranger"]], 
                                      "/raw_feature_bc_matrix"), 
                     col.names = TRUE, 
                     type = "sparse" )

# add specified metadata to object 
metadata <- read.csv(file = snakemake@params[["metadata"]], 
                     head = TRUE, 
                     sep = ",", 
                     check.names=FALSE, 
                     stringsAsFactors=FALSE, 
                     as.is=TRUE, 
                     colClasses = "character")

# if necessary, concatenate identifiers again to obtain all possible wildcards 
if(! "individual" %in% colnames(metadata)){
  metadata <- metadata %>%
    unite(col = "individual", all_of(IDENTIFIERS), sep = "_", remove = FALSE)
  print('establishing "individual" column, construct_sce_objects.R')
}

# subset data as specified by wildcard and single_cell_object_metadata_fields
metadata <- metadata |> filter(individual == individual_curr)
cols_add <- snakemake@params[["single_cell_object_metadata_fields"]]
metadata <- metadata[,colnames(metadata) %in% cols_add]

# add cell-level metadata to each barcode in the SCE object
for(i in colnames(metadata)){
  colData(sce)[i] <- rep(metadata[1,i], ncol(sce))
}

print(sce)
saveRDS(sce, file = snakemake@output[["sce_objects"]])