#-------------------------------------------------------------------------------
# in this script, results from DESeq objects are exported for every 
# combination of binary species comparisons
# for non-differentially expressed genes

library(DESeq2, quietly = TRUE)
library(tidyverse, quietly = TRUE)

#-------------------------------------------------------------------------------

tdsq_list <- base::readRDS(file = snakemake@input[["deseq_input"]])

fc_cutoff <- snakemake@params[["fc_cutoff"]]
padj_cutoff <- snakemake@params[["padj_cutoff"]]

fraction_curr <- snakemake@wildcards[["fraction"]]
species <- snakemake@params[["species"]]

# these are cell types that are so small that usually they don't work
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)

#-------------------------------------------------------------------------------
# remove objects from these cell types from the list
if(!is.null(cts_exclude)){
  tdsq_list <- tdsq_list[!names(tdsq_list) %in% cts_exclude]
}
celltypes <- names(tdsq_list)
print(celltypes)

#-------------------------------------------------------------------------------

a <- species[1]
b <- species[2]
c <- species[3]
d <- species[4]

# every combination for easier sub-setting later
comb_list<- list(
  c(a, b),
  c(a, c),
  c(a, d),
  c(b, a),
  c(b, c),
  c(b, d),
  c(c, a),
  c(c, b),
  c(c, d),
  c(d, a),
  c(d, b),
  c(d, c)
)

#-------------------------------------------------------------------------------
# RESULTS
# get results for each combination (pairwise comparison)
# specify that results should be obtained for an alternative hypothesis,
# that is, lessAbs = less than the absolute value of given log2FC threshold

res_list_list <- lapply(tdsq_list, function(tdsq){
  
  res_list <- lapply(comb_list, function(comb){
    
    res_lfcs <- DESeq2::results(tdsq,
                                lfcThreshold=fc_cutoff, 
                                altHypothesis="lessAbs", # it's log 2 FC
                                contrast=c("species", comb[1], comb[2]),
                                alpha = padj_cutoff)
    
    return(res_lfcs) 
  })
  
  for(i in 1:length(comb_list)){
    names(res_list)[[i]] <- base::paste0(comb_list[[i]][1], # a
                                         "-", 
                                         comb_list[[i]][2]) # b
  }
  
  return(res_list)
})
names(res_list_list) <- celltypes
print(names(res_list_list[[1]]))

base::saveRDS(res_list_list, file = snakemake@output[["celltype_res"]])

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# EXPORT DF
# export list of dataframes (one per cell type) for easier handling  

# a list of dataframes per cell type, each dataframe containing all comparisons
res_df_list <- lapply(res_list_list, function(res_list){
  
  res_df <- base::data.frame()
  for(j in names(res_list)){
    res_temp <- base::data.frame(res_list[[j]])
    res_temp <- tibble::rownames_to_column(res_temp, var = "gene")
    
    res_temp$comparison <- j
    res_temp$species <- stringr::str_sub(j, 1, 4)
    
    res_df <- BiocGenerics::rbind(res_df, res_temp)
  }
  
  return(res_df)
})
names(res_df_list) <- names(res_list_list)

# save list of dataframes
base::saveRDS(res_df_list, file = snakemake@output[["celltype_resdf_list"]])

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SHARED GENES
# 1. make a list of genes shared between each pairwise comp, for each cell type

res_list_shared_comp <- lapply(res_df_list, function(res_df){
  
  shared_list <- list()
  
  comparisons <- base::unique(res_df$comparison)
  for(c in comparisons){
    
    # subset by comparison
    res_df_comp <- res_df[res_df$comparison == c,]
    
    # get genes with absolute log2 FC < and padj > thresholds
    shared_genes <- res_df_comp$gene[
      which(
        base::abs(res_df_comp$log2FoldChange) < fc_cutoff & 
          res_df_comp$padj < padj_cutoff)]
    
    # export unique shared genes
    shared_genes <- base::unique(shared_genes)
    shared_list[[c]] <- shared_genes
  }
  return(shared_list)
})
print(names(res_list_shared_comp))
print(names(res_list_shared_comp[[1]]))
print(head(res_list_shared_comp[[1]][[1]]))

#-------------------------------------------------------------------------------

# 2. make a list of genes that are shared between ALL pairwise comparisons

celltype_shared_list <- lapply(res_list_shared_comp, function(shared){
  
  # for each cell type get the intersect between all pairwise comparisons
  # all appear twice (a-b and b-a), but their intersection is identical
  shared_all <- BiocGenerics::Reduce(BiocGenerics::intersect, shared)
  print(length(shared_all))
  
  return(shared_all)
})

# save list
base::saveRDS(celltype_shared_list,
              file = snakemake@output[["celltype_shared_genes_list"]])

utils::sessionInfo()
