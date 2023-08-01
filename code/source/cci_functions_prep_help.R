# helper functions for calculation of cell type interactomes

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

prepare_extraction <- function(sce, assay_use){
  
  # get the appropriate counts matrix
  if(is.null(assay_use) == TRUE){assay_use <- "logcounts"}
  
  # create countsmatrix and datasheet for extract_matrix() 
  countsmatrix <- assays(sce)[[grep(assay_use, names(assays(sce)))]]
  rownames(countsmatrix) <- rowData(sce)$ID
  
  # create datasheet
  datasheet <- data.frame(
    sample = colnames(sce),
    identity = colData(sce)$Identity,
    assignment = colData(sce)$Assignment
  )
  
  # extract identities
  idents <- levels(colData(sce)$Identity)

  return(list(
    countsmatrix = countsmatrix,
    datasheet = datasheet,
    idents = idents
  ))
}

#-------------------------------------------------------------------------------

perc_expr_genes <- function(ident, perc_df, cutoff){
  
  print(ident)
  
  expr_genes <- rownames(perc_df)[which(
    perc_df[,colnames(perc_df) == ident] >= cutoff)]
  print(length(expr_genes))
  
  return(expr_genes)
}

#-------------------------------------------------------------------------------

#perc_cells_cutoff <- function(ident_list, perc_df, cutoff){
  # list of cell types
  # pct_expr_cells = percent of expressing cells df from expr_cells_perc()
  # cutoff of minimum percentage of expressing cells required
  
#  expr_genes_pos <- which(
#    pct_expr_cells[,colnames(pct_expr_cells) == cts_list] >= cutoff)
#  return(expr_genes_pos)
#}

# returns for each cell type (= col) a vector of gene positions that are 
# expressed in a given percentage of cells per cell type as defined by cutoff

