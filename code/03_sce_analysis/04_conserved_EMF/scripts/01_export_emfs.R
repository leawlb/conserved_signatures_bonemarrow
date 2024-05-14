#-------------------------------------------------------------------------------
# export a list of marker gene dataframes and lists containing info on species-
# specificity and conservation level, one item for each cell type and fraction
# each list item contains:
# - df_cons: a dataframe containing detailed info on conserved marker genes
# - cons_emf: vector of EMFs
# - df_core: a dataframe containing detailed info on EMFs
# - ndges: vector containing all ndges of this cell type
# - marker_cons_ct: the entire dataframe containing all marker genes, conserved
#   or not, as well as information on conservation level
#TODO: remove unnecessary items after talking to Veronica
# unneccessary utems would be df_cons and df_core, replacing df_cons with a vector
# also rename all items so the naming scheme makes sense

set.seed(37)
library(tidyverse)

#-------------------------------------------------------------------------------
# load objects
load(snakemake@input[["marker_cons"]]) 
marker_cons <- markers_conservation # loaded from input

celltype_ndge_list <- readRDS(snakemake@input[["celltype_ndge_list"]])

#-------------------------------------------------------------------------------
# for each cell type, make a list item with the overlap between shared nDGEs 
# (celltype_ndge_list) and conserved marker genes (marker_cons)
cts <- names(celltype_ndge_list)
cts_list <- as.list(cts)

# make a function to collect data 
cons_emf_list <- lapply(cts_list, function(ct, 
                                           marker_cons, 
                                           celltype_ndge_list){
  
  marker_cons_ct <- as.data.frame(marker_cons[[ct]])
  shared_ct <- celltype_ndge_list[[ct]]$all
  #shared_ct <- celltype_ndge_list[[ct]]$at_least_three
  
  marker_cons_ct <- rownames_to_column(marker_cons_ct, var = "gene")
  
  marker_cons_ct$cons_marker <- vector(length = nrow(marker_cons_ct))
  marker_cons_ct$cons_marker[which(!is.na(marker_cons_ct$mmus) &
                                     !is.na(marker_cons_ct$mcas) &
                                     !is.na(marker_cons_ct$mspr) &
                                     !is.na(marker_cons_ct$mcar))] <- TRUE
  

  # only genes that were identified as markers in all four species
  marker_cons_ct_temp <- marker_cons_ct[marker_cons_ct$cons_marker == TRUE,]
  
  # get overlap
  cons_emf <- intersect(marker_cons_ct_temp$gene, shared_ct)
  print(nrow(marker_cons_ct_temp))
  print(length(cons_emf))
  
  marker_cons_ct$cons_emf <- vector(length = nrow(marker_cons_ct))
  marker_cons_ct$cons_emf[marker_cons_ct$gene %in% cons_emf] <- TRUE
  
  marker_cons_ct$ndge <- vector(length = nrow(marker_cons_ct))
  marker_cons_ct$ndge[marker_cons_ct$gene %in% shared_ct] <- TRUE
  
  # make data frames
  df_core <- marker_cons_ct_temp[marker_cons_ct_temp$gene %in% cons_emf,]
  
  df_core <- df_core[order(df_core$mmus, decreasing = TRUE),]
  marker_cons_ct_temp <- marker_cons_ct_temp[
    order(marker_cons_ct_temp$mmus, decreasing = TRUE),]

  # also add the whole thing for convenience
  return(list("df_cons" = marker_cons_ct_temp,
              "cons_emf" = cons_emf,
              "df_core" = df_core,
              "ndges" = shared_ct,
              "marker_cons_ct" = marker_cons_ct))
  
}, marker_cons,  celltype_ndge_list)
names(cons_emf_list) <- cts

saveRDS(cons_emf_list, snakemake@output[["cons_EMF_list"]])