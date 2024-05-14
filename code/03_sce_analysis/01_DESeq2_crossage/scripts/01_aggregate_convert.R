#-------------------------------------------------------------------------------

# aggregate SCE objects for each Sample (per Object_ID) and cell type, 
# then convert to Deseq2 object

library(DESeq2)
library(scuttle)
set.seed(37)

#-------------------------------------------------------------------------------
# load objects
sce <- readRDS(file = snakemake@input[["sce_input"]])

stopifnot(!is.na(sce$annotation_cluster))
stopifnot(!is.na(sce$cluster_louvain))

species <- snakemake@params["species"]

#-------------------------------------------------------------------------------

# aggregate all cells with the same combination of Object_ID and cell type 
agg <- scuttle::aggregateAcrossCells(
  sce, 
  id=colData(sce)[,c("Object_ID", "celltypes")])

cts <- levels(agg$celltypes)

#-------------------------------------------------------------------------------
# separate into list by species
agg_list <- list()
for(s in unique(agg$Species_ID)){
  agg_list[[s]] <-  agg[,agg$Species_ID == s]
}

# for each cell type, put into a separate list item to perform analysis per ct
# only keep relevant col data
dsq_list <- lapply(agg_list, function(aggl){
  dsq_list_cts <- list()
    
  for(i in cts){
    if(i %in% aggl$celltypes){
      agg <- aggl[,aggl$celltypes==i]
  
      # change colnames for compatibility and easier readability
      agg$batch <- agg$Batch_exp_day
      agg$age <- agg$Age_ID
      agg$species <- agg$Species_ID
      agg$sample <- agg$Object_ID
      agg$condition <- paste0(agg$species, "_", agg$age)

      colData(agg) <- colData(agg)[,c(
        "batch", "species", "sample", "condition",
        "age", "ncells", "Antibody_combination", "Batch_sequencing")]
    
      # convert to DeSeq2 object 
      # design will be adjusted later
      dsq <- DESeq2::DESeqDataSet(agg, design = ~ age) 
      colnames(dsq) <- agg$sample
  
      # factor for nicer visualisations
      colData(dsq)$age <- factor(colData(dsq)$age, levels = c("yng", "old"))
      colData(dsq)$Antibody_combination <- factor(
        colData(dsq)$Antibody_combination,
        levels = unique(colData(dsq)$Antibody_combination))
      colData(dsq)$sample <- factor(colData(dsq)$sample,
                                    levels = unique(colData(dsq)$sample))
    
      dsq_list_cts[[i]] <- dsq
    }else{
      dsq_list_cts[[i]] <- NULL
    }
  }
  return(dsq_list_cts)
})
# each list item contains an object for each cell type, containing all 
# samples from both ages per species

#-------------------------------------------------------------------------------

output_paths <- snakemake@output[["deseq_output"]]
print(output_paths)

# save by species list
saveRDS(dsq_list$mmus, file = output_paths[[grep("mmus", output_paths)]])
saveRDS(dsq_list$mcas, file = output_paths[[grep("mcas", output_paths)]])
saveRDS(dsq_list$mcar, file = output_paths[[grep("mcar", output_paths)]])
saveRDS(dsq_list$mspr, file = output_paths[[grep("mspr", output_paths)]])

sessionInfo()