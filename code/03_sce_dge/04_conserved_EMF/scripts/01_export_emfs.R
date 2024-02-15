#-------------------------------------------------------------------------------
set.seed(37)
library(tidyverse)

# load objects
load(snakemake@input[["marker_cons"]])
marker_cons <- markers_conservation # loaded from input

celltype_res_list_shared <- readRDS(snakemake@input[["celltype_res_list_shared"]])
#celltype_res_list_shared <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/05_nres/PC_0.05_FC_1.5/res_hsc_celltype_shared")

cts <- names(celltype_res_list_shared)
cts2 <- names(marker_cons)

# same cell type names but different sequence
table(is.na(match(cts, cts2)))
table(is.na(match(cts2, cts)))

# for each cell type, make a list item with the overlapt between shared nDGEs 
# (celltype_res_list_shared) and conserved marker genes (marker_cons)
cts_list <- as.list(cts)

cons_emf_list <- lapply(cts_list, function(ct, marker_cons, 
                                            celltype_res_list_shared){
  
  marker_cons_ct <- as.data.frame(marker_cons[[ct]])
  shared_ct <- celltype_res_list_shared[[ct]]$all
  #shared_ct <- celltype_res_list_shared[[ct]]$at_least_three
  
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
  
}, marker_cons,  celltype_res_list_shared)
names(cons_emf_list) <- cts

saveRDS(cons_emf_list, snakemake@output[["cons_EMF_list"]])