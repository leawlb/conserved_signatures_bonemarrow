#-------------------------------------------------------------------------------
# aggregate SCE objects for each Sample, then convert to Deseq2 object

library(DESeq2, quietly = TRUE)
library(scuttle, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# load objects
sce <- base::readRDS(file = snakemake@input[["sce_input"]])

stopifnot(!is.na(sce$annotation_cluster))
stopifnot(!is.na(sce$cluster_louvain))

#-------------------------------------------------------------------------------
# aggregate all cells with the same combination of Object_ID (sample) and
# annotation_cluster (preliminary cell type)
# using unaltered counts before batch correction
agg <- scuttle::aggregateAcrossCells(
  sce,
  use.assay.type = "counts_before_BC",
  id = colData(sce)[,c("Object_ID", "annotation_cluster")])

print(agg)
print(levels(agg$annotation_cluster))

#-------------------------------------------------------------------------------
# put each cell type object into a separate list item to perform analysis per ct
# only keep relevant col data

agg_list <- list()
for(i in levels(agg$annotation_cluster)){
  agg_list[[i]] <- agg[,agg$annotation_cluster==i]
  
  agg_list[[i]]$batch <- factor(agg_list[[i]]$Batch_exp_day)
  agg_list[[i]]$condition <- factor(agg_list[[i]]$Species_ID) # SET CONDITION
  agg_list[[i]]$sample <- factor(agg_list[[i]]$Object_ID)
  agg_list[[i]]$age <- factor(agg_list[[i]]$Age_ID)
  
  colData(agg_list[[i]]) <- colData(agg_list[[i]])[,c("batch", 
                                                      "condition", 
                                                      "sample",
                                                      "age",
                                                      "annotation_cluster",
                                                      "ncells", 
                                                      "Antibody_combination", 
                                                      "Batch_sequencing")]
  
  print(i)
  print(agg_list[[i]])
}
# each list item contains an object for each cell type, containing all 
# samples from each species and age

#-------------------------------------------------------------------------------
# convert to DESeq2 object (sequence in design will be adjusted in 07)
# "batch" and "age" are covariates, "condition" the condition to be tested 
# = species

dsq_list <- lapply(agg_list, function(agg){
  
  print(agg)
  print(agg$annotation_cluster[1])
  dsq <- DESeq2::DESeqDataSet(agg, design = ~ batch + age + condition) 
  colnames(dsq) <- agg$sample
  
  # factor columns
  colData(dsq)$condition <- factor(
    colData(dsq)$condition, 
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
