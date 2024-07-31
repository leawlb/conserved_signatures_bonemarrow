#-------------------------------------------------------------------------------
# aggregate SCE objects for each Sample (per Object_ID) and cell type, 
# then convert to Deseq2 object

library(DESeq2, quietly = TRUE)
library(scuttle, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# load objects
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

stopifnot(!is.na(sce$annotation_cluster))
stopifnot(!is.na(sce$cluster_louvain))
stopifnot(!is.na(sce$celltypes))

species <- snakemake@params["species"]
ct_exclude <- snakemake@params["cell_types_exclude"]
print(ct_exclude)

sce <- sce[,!sce$celltypes %in% ct_exclude]
sce$celltypes <- factor(
  sce$celltype, 
  levels = levels(sce$celltypes)[!levels(sce$celltypes) %in% ct_exclude])
print(levels(sce$celltypes))

#-------------------------------------------------------------------------------

# aggregate all cells with the same combination of Object_ID and cell type 
agg <- scuttle::aggregateAcrossCells(
  sce, 
  use.assay.type = "counts_before_BC", # raw counts before batch correction
  id=colData(sce)[,c("Object_ID", "celltypes")])

cts <- levels(agg$celltypes)

#-------------------------------------------------------------------------------
# separate into list by species
agg_list <- list()
for(s in base::unique(agg$Species_ID)){
  agg_list[[s]] <- agg[,agg$Species_ID == s]
}

# for each cell type, put into a separate list item to perform analysis per ct
# only keep relevant col data
dsq_list <- lapply(agg_list, function(aggl){
  dsq_list_cts <- list()
    
  for(ct in cts){
    if(ct %in% aggl$celltypes){
      
      agg <- aggl[,aggl$celltypes==ct]
  
      # change colnames for compatibility and easier readability
      agg$batch <- agg$Batch_exp_day
      agg$age <- agg$Age_ID
      agg$species <- agg$Species_ID
      agg$sample <- agg$Object_ID
      agg$condition <- base::paste0(agg$species, "_", agg$age)

      colData(agg) <- colData(agg)[,c(
        "batch", 
        "species", 
        "sample", 
        "condition",
        "age",
        "ncells",
        "Antibody_combination",
        "Batch_sequencing")]
    
      # convert to DeSeq2 object 
      # design will be adjusted later
      print(ct)
      dsq <- DESeq2::DESeqDataSet(agg, design = ~ age) 
      colnames(dsq) <- agg$sample
  
      # factor for nicer visualisation
      colData(dsq)$age <- factor(
        colData(dsq)$age, 
        levels = base::unique(colData(dsq)$age))
      colData(dsq)$Antibody_combination <- factor(
        colData(dsq)$Antibody_combination,
        levels = base::unique(colData(dsq)$Antibody_combination))
      colData(dsq)$sample <- factor(
        colData(dsq)$sample,
        levels = base::unique(colData(dsq)$sample))
    
      dsq_list_cts[[ct]] <- dsq
    }else{
      dsq_list_cts[[ct]] <- NULL
    }
  }
  return(dsq_list_cts)
})
# each list item (per species) contains an object for each cell type, 
# containing all samples from both ages 

#-------------------------------------------------------------------------------

output_paths <- snakemake@output[["deseq_output"]]
print(output_paths)

# save by species list
base::saveRDS(dsq_list$mmus, 
              file = output_paths[[base::grep("mmus", output_paths)]])
base::saveRDS(dsq_list$mcas,
              file = output_paths[[base::grep("mcas", output_paths)]])
base::saveRDS(dsq_list$mcar, 
              file = output_paths[[base::grep("mcar", output_paths)]])
base::saveRDS(dsq_list$mspr, 
              file = output_paths[[base::grep("mspr", output_paths)]])

utils::sessionInfo()