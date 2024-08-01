#-------------------------------------------------------------------------------
# aggregate SCE objects for each Sample, then convert to Deseq2 object

library(DESeq2, quietly = TRUE)
library(scuttle, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- base::readRDS(file = snakemake@input[["sce_input"]])

stopifnot(!is.na(sce$annotation_cluster))
stopifnot(!is.na(sce$cluster_louvain))

#-------------------------------------------------------------------------------
# aggregate all cells with the same combination of Object_ID (sample) and
# cell type
# using unaltered counts before batch correction
agg <- scuttle::aggregateAcrossCells(
  sce, 
  use.assay.type = "counts_before_BC",
  id=colData(sce)[,c("Object_ID", "celltypes")])

print(agg)

#-------------------------------------------------------------------------------
# for each cell type, put into a separate list item to perform analysis per ct
# only keep relevant col data

agg_list <- list()
for(ct in levels(agg$celltypes)){
  
  print(ct)
  agg_list[[ct]] <- agg[,agg$celltypes==ct]
  
  agg_list[[ct]]$batch <- factor(agg_list[[ct]]$Batch_exp_day)
  agg_list[[ct]]$species <- factor(agg_list[[ct]]$Species_ID)
  agg_list[[ct]]$sample <- factor(agg_list[[ct]]$Object_ID)
  agg_list[[ct]]$age <- factor(agg_list[[ct]]$Age_ID)
  
  colData(agg_list[[ct]]) <- colData(agg_list[[ct]])[,c("batch",
                                                        "species", 
                                                        "sample",
                                                        "age", 
                                                        "ncells", 
                                                        "Antibody_combination", 
                                                        "Batch_sequencing")]
}
# each list item contains an object for each cell type, containing all 
# samples from each species and age

#-------------------------------------------------------------------------------
# convert to DESeq2 object (sequence in design will be adjusted in 07)
# "batch" and "age" are covariates, "condition" the condition to be tested 
# = species

dsq_list <- lapply(agg_list, function(agg){
  
  print(agg)
  
  dsq <- DESeq2::DESeqDataSet(agg, design = ~ batch + age + species) 
  colnames(dsq) <- agg$sample
  
  # factor columns
  colData(dsq)$species <- factor(
    colData(dsq)$species,
    levels = c("mmus", "mcas", "mspr", "mcar"))
  colData(dsq)$age <- factor(
    colData(dsq)$age,
    levels = c("yng", "old"))
  colData(dsq)$Antibody_combination <- factor(
    colData(dsq)$Antibody_combination,
    levels = base::unique(colData(dsq)$Antibody_combination))
  colData(dsq)$sample <- factor(
    colData(dsq)$sample,
    levels = base::unique(colData(dsq)$sample))
  
  return(dsq)
})

#-------------------------------------------------------------------------------
base::saveRDS(dsq_list, file = snakemake@output[["deseq_output"]])

utils::sessionInfo()