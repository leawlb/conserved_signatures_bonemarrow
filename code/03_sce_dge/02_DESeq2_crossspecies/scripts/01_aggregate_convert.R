
# aggregate SCE objects for each Sample, then convert to Deseq2 object

library(DESeq2)
library(scuttle)
library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------
# load objects
sce <- readRDS(file = snakemake@input[["sce_input"]])

stopifnot(!is.na(sce$annotation_cluster))
stopifnot(!is.na(sce$cluster_louvain))

#-------------------------------------------------------------------------------
# aggregate all cells with the same combination of Object ID and cell type 
agg <- aggregateAcrossCells(sce, id=colData(sce)[,c("Object_ID",
                                                    "celltypes")])

print(agg)
colData(agg)

# for each cell type, put into a separate list item to perform analysis per ct
# only keep relevant col data
agg_list <- list()
for(i in levels(agg$celltypes)){
  agg_list[[i]] <- agg[,agg$celltypes==i]
  
  agg_list[[i]]$batch <- factor(agg_list[[i]]$Batch_exp_day)
  agg_list[[i]]$species <- factor(agg_list[[i]]$Species_ID)
  agg_list[[i]]$sample <- factor(agg_list[[i]]$Object_ID)
  agg_list[[i]]$age <- factor(agg_list[[i]]$Age_ID)
  
  colData(agg_list[[i]]) <- colData(agg_list[[i]])[,c("batch", "species", 
                                                      "sample", "age", "ncells", 
                                                      "Antibody_combination", 
                                                      "Batch_sequencing")]
  
}
# each list item contains an object for each cell type, containing all 
# samples from each species and age

#-------------------------------------------------------------------------------
# convert to DeSeq2 object (sequence in design is not too important for now)
# "batch" and "age" are covariates, "species" the species to be tested

dsq_list <- lapply(agg_list, function(agg){
  dsq <- DESeqDataSet(agg, design = ~ batch + age + species) 
  colnames(dsq) <- agg$sample
  
  colData(dsq)$species <- factor(colData(dsq)$species,
                                   levels = c("mmus", "mcas", "mspr", "mcar"))
  colData(dsq)$age <- factor(colData(dsq)$age, levels = c("yng", "old"))
  colData(dsq)$Antibody_combination <- factor(
    colData(dsq)$Antibody_combination,
    levels = unique(colData(dsq)$Antibody_combination))
  colData(dsq)$sample <- factor(colData(dsq)$sample,
                                levels = unique(colData(dsq)$sample))
  
  return(dsq)
})

saveRDS(dsq_list, file = snakemake@output[["deseq_output"]])