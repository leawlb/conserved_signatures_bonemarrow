
library(DESeq2)
library(scuttle)
library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------
# load objects
sce <- readRDS(file = snakemake@input[["sce_11"]])

#-------------------------------------------------------------------------------
# aggregate all cells with the same combination of Object ID and cell type 
agg <- aggregateAcrossCells(sce, id=colData(sce)[,c("Object_ID",
                                                    "cluster_louvain")])

print(agg)
# for each cell type, put into a separate list item to perform analysis per ct
# only keep relevant col data
agg_list <- list()
for(i in levels(agg$cluster_louvain)){
  agg_list[[i]] <- agg[,agg$cluster_louvain==i]
  
  agg_list[[i]]$batch <- factor(agg_list[[i]]$Batch_exp_day)
  agg_list[[i]]$condition <- factor(agg_list[[i]]$Species_ID)
  agg_list[[i]]$sample <- factor(agg_list[[i]]$Object_ID)
  agg_list[[i]]$age <- factor(agg_list[[i]]$Age_ID)
  
  colData(agg_list[[i]]) <- colData(agg_list[[i]])[,c("batch", "condition", 
                                                      "sample", "age", "ncells", 
                                                      "Antibody_combination", 
                                                      "Batch_sequencing")]
  
}

#-------------------------------------------------------------------------------
# convert to DeSeq2 object (sequence in design is not too important for now)
# "batch" and "age" are covariates, "condition" the condition to be tested

dsq_list <- lapply(agg_list, function(agg){
  dsq <- DESeqDataSet(agg, design = ~ batch + age + condition) 
  colnames(dsq) <- agg$sample
  
  colData(dsq)$condition <- factor(colData(dsq)$condition,
                                   levels = c("mmus", "mcas", "mspr", "mcar"))
  colData(dsq)$age <- factor(colData(dsq)$age, levels = c("yng", "old"))
  colData(dsq)$Antibody_combination <- factor(colData(dsq)$Antibody_combination,
                                              levels = unique(colData(dsq)$Antibody_combination))
  colData(dsq)$sample <- factor(colData(dsq)$sample,
                                levels = unique(colData(dsq)$sample))
  
  return(dsq)
})

saveRDS(dsq_list, file = snakemake@output[["dsq_11"]])
