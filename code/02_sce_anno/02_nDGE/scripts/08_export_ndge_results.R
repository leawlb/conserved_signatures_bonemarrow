#-------------------------------------------------------------------------------
# results from DESeq objects (nDGES) are exported for every cluster

library(DESeq2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------

tdsq_list <- base::readRDS(file = snakemake@input[["deseq_input"]])

fc_cutoff <- snakemake@params[["fc_cutoff"]]
padj_cutoff <- snakemake@params[["padj_cutoff"]]

fraction_curr <- snakemake@wildcards[["fraction"]]
species <- snakemake@params[["species"]]

clusters <- names(tdsq_list)
print(clusters)

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
                                contrast=c("condition", comb[1], comb[2]),
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
names(res_list_list) <- clusters
names(res_list_list[[1]])

base::saveRDS(res_list_list, file = snakemake@output[["cluster_res"]])

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
base::saveRDS(res_df_list, file = snakemake@output[["cluster_resdf_list"]])

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
    
    # get genes with absolute log2 FC < and padj < thresholds
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

res_list_shared_all <- lapply(res_list_shared_comp, function(shared){
  
  # for each cell type get the intersect between all pairwise comparisons
  # all appear twice (a-b and b-a), but their intersection is identical
  shared_all <- BiocGenerics::Reduce(BiocGenerics::intersect, shared)
  
  return(shared_all)
})

#-------------------------------------------------------------------------------

# 3. make a list of genes that are shared between three of four species

res_list_shared_three <- lapply(res_list_shared_comp, function(shared){

  print(names(shared))
  
  # manually get intersections between required pairwise comparisons
  # a b c
  shared_between_mmus_mcas_mspr <- BiocGenerics::Reduce(
    BiocGenerics::intersect, 
    list(shared$`mmus-mcas`, # a b 
         shared$`mmus-mspr`, # a c
         shared$`mcas-mmus`, # b a 
         shared$`mcas-mspr`, # b c
         shared$`mspr-mmus`, # c a 
         shared$`mspr-mcas`  # c b = a b c 
  ))
  
  # a b d
  shared_between_mmus_mcas_mcar <- BiocGenerics::Reduce(
    BiocGenerics::intersect, 
    list(shared$`mmus-mcas`,
         shared$`mmus-mcar`,
         shared$`mcas-mmus`,
         shared$`mcas-mcar`,
         shared$`mcar-mmus`,
         shared$`mcar-mcas`
    ))
  
  # a c d
  shared_between_mmus_mspr_mcar <- BiocGenerics::Reduce(
    BiocGenerics::intersect, 
    list(shared$`mmus-mspr`,
         shared$`mmus-mcar`,
         shared$`mspr-mmus`,
         shared$`mspr-mcar`,
         shared$`mcar-mmus`,
         shared$`mcar-mspr`
    ))
  
  # b c d
  shared_between_mcas_mspr_mcar <- BiocGenerics::Reduce(
    BiocGenerics::intersect, 
    list(shared$`mcas-mspr`,
         shared$`mcas-mcar`,
         shared$`mspr-mcas`,
         shared$`mspr-mcar`,
         shared$`mcar-mcas`,
         shared$`mcar-mspr`
    ))  
  
  # do some random tests
  stopifnot(shared$`mcas-mspr` == shared$`mspr-mcas`)
  stopifnot(shared$`mmus-mspr` == shared$`mspr-mmus`)
  stopifnot(shared$`mcar-mcas` == shared$`mcas-mcar`)
  
  # add all genes together, as all these genes have already passed the
  # requirements
  genes_three <- base::unique(c(
    shared_between_mmus_mcas_mspr,
    shared_between_mmus_mcas_mcar,
    shared_between_mmus_mspr_mcar,
    shared_between_mcas_mspr_mcar
    
  ))  
  
  print(length(genes_three))
  return(genes_three)
  
})

export_list <- list()
for(ct in names(res_list_shared_three)){
  export_list[[ct]][["all"]] <- res_list_shared_all[[ct]]
  export_list[[ct]][["three"]] <- res_list_shared_three[[ct]]
}

# save list
base::saveRDS(export_list, file = snakemake@output[["cluster_sharedgenes_list"]])

utils::sessionInfo()
