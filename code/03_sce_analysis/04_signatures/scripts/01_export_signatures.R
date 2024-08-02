#-------------------------------------------------------------------------------
# export a list of marker gene data frames and lists containing info on 
# - conserved marker genes
# - conserved signature genes
# - nDGEs
# - genes used for subclustering
# one item for each cell type and fraction

set.seed(37)
library(tidyverse)

#-------------------------------------------------------------------------------
# load objects
# conserved markers (from Veronica)
load(snakemake@input[["marker_cons"]]) 
marker_cons <- markers_conservation # loaded from input

# non-differentially expressed genes
ct_ndge_list <- base::readRDS(snakemake@input[["celltype_ndge_list"]])

# genes used for subclustering from some cell types
# these might be excluded from downstream analysis as they were used 
# for semi-supervised clustering/cell type identification
fraction_curr <- snakemake@wildcards[["fraction"]]

subclustering_genes_path <- snakemake@input[["subclustering_genes"]]
print(subclustering_genes_path)
genes_subcl_df <- utils::read.csv(file = subclustering_genes_path, 
                                  header = TRUE, 
                                  sep = ";", 
                                  check.names=FALSE, 
                                  stringsAsFactors=FALSE, 
                                  as.is=TRUE, 
                                  colClasses = "character")
print(head(genes_subcl_df))
genes_subcl_df <- genes_subcl_df[genes_subcl_df$fraction == fraction_curr,]
genes_subcl_df <- genes_subcl_df[genes_subcl_df$purpose == "subclustering",]
genes_subcl <- genes_subcl_df$gene

#-------------------------------------------------------------------------------

cts <- names(ct_ndge_list)
cts <- cts[cts %in% names(marker_cons)]
cts_list <- as.list(cts)



# make a function to collect data 
cons_signature_list <- lapply(cts_list, function(ct){
  
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
names(cons_signature_list) <- cts

base::saveRDS(cons_signature_list, snakemake@output[["signature_list"]])

utils::sessionInfo()