#-------------------------------------------------------------------------------
# export a list of marker gene data frames and lists containing info on 
# - conserved signature genes
# - conserved marker genes
# - nDGEs
# - all BL6 marker genes
# - genes used for subclustering
# one item for each cell type and fraction

set.seed(37)

library(tidyverse, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# load

fraction_curr <- snakemake@wildcards[["fraction"]]

# conserved markers 
load(snakemake@input[["marker_cons"]]) 
if(fraction_curr == "hsc"){
  marker_cons <- markers_conservation_hsc # loaded from input
}else if(fraction_curr == "str"){
  marker_cons <- markers_conservation_str # loaded from input
}

# non-differentially expressed genes
ct_ndge_list <- base::readRDS(snakemake@input[["celltype_ndge_list"]])

# cell types that were not used to extract signature are excluded because
# they cannot be separated after excluding them
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)

#-------------------------------------------------------------------------------

# sce contains info on subclustering genes
# genes used for sub-clustering from some cell types might be excluded from 
# downstream analysis as they were used to define cell types
sce <- base::readRDS(snakemake@input[["sce_input"]])
print(sce)

#-------------------------------------------------------------------------------

cts <- names(ct_ndge_list)
cts <- cts[cts %in% names(marker_cons)]
cts <- cts[!cts %in% cts_exclude]
cts_list <- as.list(cts)

# make a function to collect data 
geneset_list <- lapply(cts_list, function(ct){
  
  marker_cons_ct <- base::as.data.frame(marker_cons[[ct]])
  shared_ct <- ct_ndge_list[[ct]]

  # prepare df
  marker_cons_ct <- tibble::rownames_to_column(marker_cons_ct, var = "gene")
  
  # get all genes that are conserved markers in all four species
  marker_cons_ct$conserved_marker <- vector(length = nrow(marker_cons_ct))
  marker_cons_ct$conserved_marker[which(!is.na(marker_cons_ct$mmus) &                                         
                                          !is.na(marker_cons_ct$mcas) &
                                          !is.na(marker_cons_ct$mspr) &
                                          !is.na(marker_cons_ct$mcar))] <- TRUE
  
  # get conserved signature = overlap between conserved markers and nDGEs
  marker_cons_ct_temp <- marker_cons_ct[
    marker_cons_ct$conserved_marker == TRUE,]
  conserved_signature <- BiocGenerics::intersect(marker_cons_ct_temp$gene, 
                                                 shared_ct)
  # add info on non-differentially expressed genes to df
  marker_cons_ct$ndge <- vector(length = nrow(marker_cons_ct))
  marker_cons_ct$ndge[
    marker_cons_ct$gene %in% shared_ct] <- TRUE
  
  # add info on conserved signature to df
  marker_cons_ct$conserved_signature <- vector(length = nrow(marker_cons_ct))
  marker_cons_ct$conserved_signature[
    marker_cons_ct$gene %in% conserved_signature] <- TRUE

  # add info on genes used for subclustering to df
  # cell type names are chosen so that they cannot grep another cell type
  print(ct)
  genes_subcl <- rowData(sce)$Symbol[
    base::grep(ct, rowData(sce)$subclustering_genes)]
  if(length(genes_subcl) == 0){
    genes_subcl <- NULL
  }
  print(genes_subcl)
  
  marker_cons_ct$genes_subclustering <- vector(length = nrow(marker_cons_ct))
  marker_cons_ct$genes_subclustering[
    marker_cons_ct$gene %in% genes_subcl] <- TRUE
  
  # also add the whole thing for convenience
  return(list("conserved_df" = marker_cons_ct,
              "conserved_markers" = marker_cons_ct_temp$gene,
              "conserved_signature" = conserved_signature,
              "ndges" = shared_ct,
              "genes_subclustering" = genes_subcl
              ))
  
})
names(geneset_list) <- cts

#-------------------------------------------------------------------------------

base::saveRDS(geneset_list, snakemake@output[["geneset_list"]])

utils::sessionInfo()
