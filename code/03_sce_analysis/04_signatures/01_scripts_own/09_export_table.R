
#-------------------------------------------------------------------------------
# export signature genes as a useful table
# to be used as supplementary table for the manuscript and as input for
# other analyses

#-------------------------------------------------------------------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(37)

#-------------------------------------------------------------------------------

library(dplyr)

#-------------------------------------------------------------------------------
# load

geneset_list <- base::readRDS(snakemake@input[["geneset_list"]])

#-------------------------------------------------------------------------------

# add the cell type to each dataframe
for(ct in names(geneset_list)){
  geneset_list[[ct]]$ct <- ct
}

# generate a dataframe for export
make_export_df <- function(list){
  
  genes_sub <- list$genes_subclustering
  genes_sign <- list$conserved_signature
  
  temp_df <- data.frame(
    "signature_gene" = genes_sign,
    "celltype" =  base::rep(list$ct, length(genes_sign)),
    "subclustering_gene_for_ct" = vector(length = length(genes_sign))
  )
  
  # add info on whether a gene was a sub-clustering gene for the cell type
  temp_df$subclustering_gene_for_ct[which(temp_df$signature_gene %in% genes_sub)] <- TRUE
  return(temp_df)
}

export_df_list <- lapply(geneset_list, make_export_df)
export_df <- dplyr::bind_rows(export_df_list)
export_df <- export_df[order(export_df$signature_gene, decreasing = FALSE),]

head(export_df)

#-------------------------------------------------------------------------------

readr::write_delim(
  export_df, 
  snakemake@output[["signature_table"]], 
  delim = ";", 
  append = FALSE)
