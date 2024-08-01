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
  
  return(shared_all)
})

# save list
base::saveRDS(celltype_shared_list,
              file = snakemake@output[["celltype_shared_genes_list"]])

utils::sessionInfo()

# #-------------------------------------------------------------------------------
# # RESULTS
# # get results for each combination
# # specify that results should be obtained for an alternative hypothesis,
# # that is, less Abs = less than the absolute value of logFC threshold
# 
# res_list_list <- lapply(tdsq_list, function(tdsq){
#   
#   res_list <- lapply(comb_list, function(comb){
#     res_lfcs <- DESeq2::results(tdsq,
#                                 lfcThreshold=fc_cutoff, altHypothesis="lessAbs",
#                                 contrast=c("species", comb[1], comb[2]),
#                                 alpha = padj_cutoff)
#     
#     return(res_lfcs)
#   })
#   for(i in 1:length(comb_list)){
#     names(res_list)[[i]] <- paste0(comb_list[[i]][1], "-", comb_list[[i]][2])
#   }
#   return(res_list)
# })
# names(res_list_list) <- celltypes
# names(res_list_list[[1]])
# 
# saveRDS(res_list_list, file = snakemake@output[["celltype_res"]])
# 
# #-------------------------------------------------------------------------------
# # Get results DF
# 
# # add info
# res_df_list <- lapply(res_list_list, function(res_list){
#   
#   res_df <- data.frame()
#   for(j in names(res_list)){
#     res_temp <- data.frame(res_list[[j]])
#     res_temp <- rownames_to_column(res_temp, var = "gene")
# 
#     res_temp$comparison <- j
#     res_temp$species <- str_sub(j, 1, 4)
# 
#     res_df <- rbind(res_df, res_temp)
#   }
#   return(res_df)
# })
# names(res_df_list) <- names(res_list_list)
# saveRDS(res_df_list, file = snakemake@output[["celltype_res_dfs"]])
# 
# #-------------------------------------------------------------------------------
# # make a list of genes shared between each species for each cell type
# 
# # these are genes that are shared in the corresponding comparisons
# res_list_shared_comp <- lapply(res_df_list, function(res_df){
# 
#   shared_list <- list()
#   for(i in unique(res_df$comparison)){
#     res_df_temp <- res_df[res_df$comparison == i,]
#     shared_genes <- res_df_temp$gene[
#       which(abs(res_df_temp$log2FoldChange) < fc_cutoff & 
#               res_df_temp$padj < padj_cutoff)]
#     shared_genes <- unique(shared_genes)
#     shared_list[[i]] <- unique(shared_genes)
#   }
#   
#   return(shared_list)
# })
# 
# # these are genes that are shared between all, or three of four species
# res_list_shared_fin <- lapply(res_list_shared_comp, function(shared){
#   
#   return_list <- list()
#   
#   for(i in species){
#     print(i)
#     comp_names <- names(shared)[grep(i, str_sub(names(shared), 1, 4))]
#     print(comp_names)
#     
#     genes_all <- unique(c(shared[[comp_names[1]]], 
#                           shared[[comp_names[2]]], 
#                           shared[[comp_names[3]]]))
#     #length(genes_all)
#     
#     #length(which(genes_all %in% shared[[comp_names[1]]]))
#     #length(which(genes_all %in% shared[[comp_names[2]]]))
#     #length(which(genes_all %in% shared[[comp_names[3]]]))
#     
#     # these are genes that are shared between species i and one of the three 
#     # other species
#     shared_12 <- intersect(shared[[1]], shared[[2]])
#     shared_13 <- intersect(shared[[1]], shared[[3]])
#     shared_32 <- intersect(shared[[3]], shared[[2]])
#     #length(shared_12)
#     #length(shared_13)
#     #length(shared_32)
#     
#     # these are genes that are shared between species i and all others
#     shared_all <- intersect(shared_12, shared[[3]])
#     #length(shared_all)
#     
#     return_list[[i]][[paste0(comp_names[1], "_", comp_names[2])]] <- shared_12
#     return_list[[i]][[paste0(comp_names[1], "_", comp_names[3])]] <- shared_13
#     return_list[[i]][[paste0(comp_names[3], "_", comp_names[2])]] <- shared_32
#     return_list[[i]][["all"]] <- shared_all
#   }
#   return(return_list)
# })
# 
# # extract only the intersection between all these genes for export
# genes_three_species_list <- lapply(res_list_shared_fin, function(res_list){
# 
#   # this will contain many duplicates 
#   # because e.g. mmus_mspr and mspr_mmus are redundant
#   genes_three <- unique(c(
#     # all genes shared between mmus and two more species
#     res_list$mmus$`mmus-mspr_mmus-mcas`, res_list$mmus$`mmus-mspr_mmus-mcar`,
#     res_list$mmus$`mmus-mcar_mmus-mcas`,
#     # all genes shared between mcas and two more species
#     res_list$mcas$`mcas-mmus_mcas-mspr`, res_list$mcas$`mcas-mmus_mcas-mcar`,
#     res_list$mcas$`mcas-mcar_mcas-mspr`,
#     # all genes shared between mspr and two more species
#     res_list$mspr$`mspr-mmus_mspr-mcas`, res_list$mspr$`mspr-mmus_mspr-mcar`,
#     res_list$mspr$`mspr-mcar_mspr-mcas`,
#     # all genes shared between mcar and two more species
#     res_list$mcar$`mcar-mmus_mcar-mspr`, res_list$mcar$`mcar-mmus_mcar-mcas`,
#     res_list$mcar$`mcar-mcas_mcar-mspr`
#   
#   ))  
#   
#   print("checking quality")
#   print(length(genes_three))
#   
#   #print(identical(res_list$mmus$all, res_list$mcas$all))  
#   #print(identical(res_list$mmus$all, res_list$mspr$all))  
#   #print(identical(res_list$mmus$all, res_list$mcar$all))  
#   
#   # all shared genes are in genes_three but not the other way around
#   #print(table(genes_three %in% res_list$mmus$all))
#   #print(table(res_list$mmus$all %in% genes_three))
#   
#   return_list <- list(res_list$mmus$all, genes_three)
#   names(return_list) <- c("all", "at_least_three")
#   return(return_list)
# 
# })
# 
# saveRDS(genes_three_species_list, 
#         file = snakemake@output[["celltype_res_list_shared"]])
# 
# #-------------------------------------------------------------------------------
# # REARRANGE 
# # rearrange by species for species marker csv export and species-specificity
# 
# species_res_list_list <- list()
# for(s in species){
#   for(c in celltypes){
#     
#     species_res_list_list[[s]][[c]] <- res_list_list[[c]][
#       grep(s, names(res_list_list[[c]]))]
#   }
# }
# names(species_res_list_list) 
# names(species_res_list_list[[1]])
# 
# saveRDS(species_res_list_list, file = snakemake@output[["species_res"]])
