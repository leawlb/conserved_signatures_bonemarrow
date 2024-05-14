#-------------------------------------------------------------------------------
# aggregate SCE objects for each Sample, then convert to Deseq2 object

library(DESeq2, quietly = TRUE)
library(scuttle, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------
# load objects
sce <- readRDS(file = snakemake@input[["sce_input"]])

stopifnot(!is.na(sce$annotation_cluster))
stopifnot(!is.na(sce$cluster_louvain))

#-------------------------------------------------------------------------------
# aggregate all cells with the same combination of Object ID and cell type 
agg <- scuttle::aggregateAcrossCells(sce,
                                     id=colData(sce)[,c("Object_ID",
                                                        "annotation_cluster")])
print(agg)

# for each cell type, put into a separate list item to perform analysis per ct
# only keep relevant col data
agg_list <- list()
for(i in levels(agg$annotation_cluster)){
  agg_list[[i]] <- agg[,agg$annotation_cluster==i]
  
  agg_list[[i]]$batch <- factor(agg_list[[i]]$Batch_exp_day)
  agg_list[[i]]$condition <- factor(agg_list[[i]]$Species_ID)
  agg_list[[i]]$sample <- factor(agg_list[[i]]$Object_ID)
  agg_list[[i]]$age <- factor(agg_list[[i]]$Age_ID)
  
  colData(agg_list[[i]]) <- colData(agg_list[[i]])[,c("batch", "condition", 
                                                      "sample", "age", "ncells", 
                                                      "Antibody_combination", 
                                                      "Batch_sequencing")]
}
# each list item contains an object for each cell type, containing all 
# samples from each species and age

#-------------------------------------------------------------------------------
# convert to DeSeq2 object (sequence in design is not too important for now)
# "batch" and "age" are covariates, "condition" the condition to be tested

dsq_list <- lapply(agg_list, function(agg){
  dsq <- DESeq2::DESeqDataSet(agg, design = ~ batch + age + condition) 
  colnames(dsq) <- agg$sample
  
  colData(dsq)$condition <- factor(colData(dsq)$condition,
                                   levels = c("mmus", "mcas", "mspr", "mcar"))
  colData(dsq)$age <- factor(colData(dsq)$age, levels = c("yng", "old"))
  colData(dsq)$Antibody_combination <- factor(
    colData(dsq)$Antibody_combination,
    levels = unique(colData(dsq)$Antibody_combination))
  colData(dsq)$sample <- factor(colData(dsq)$sample,
                                levels = unique(colData(dsq)$sample))
  
  return(dsq)
})

#-------------------------------------------------------------------------------
saveRDS(dsq_list, file = snakemake@output[["deseq_output"]])

sessionInfo()
