# functions for calculation of cell type interactomes

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------

## PERCENTAGE OF EXPRESSED CELLS (expr_cells_perc)

#TODO: APPROVED!

expr_cells_perc <- function(sce, assay_use = assay_use){

  pct_expr_cells_df <- data.frame(
    row.names = rowData(sce)$ENSMUS_ID)

  idents <- levels(sce$Identity)
  ident_list <- as.list(idents)
  perc_list <- lapply(ident_list, function(x){
    
    perc <- vector()
    # subset sce objects by identity 
    sce_x <- sce[,colData(sce)$Identity == x]

    print(paste("........ ", x, ": ", ncol(sce_x), sep = "")) #track progress
    # for each gene i
    for(i in 1:nrow(pct_expr_cells_df)){
      # get all cells of identity x that express gene i
      nr_cell_expr <- length(which(assays(sce_x)[[
        grep(assay_use, names(assays(sce_x)))]][i,] > 0))
      
      # divide by total number of cells per identity to get percentage
      perc[i] <- nr_cell_expr/ncol(sce_x)
    }
    return(perc)
  })
  
  for(i in 1:length(perc_list)){
    pct_expr_cells_df[,i] <- perc_list[[i]]
    colnames(pct_expr_cells_df)[i] <- idents[i]
  }
  pct_expr_cells_df$gene_symbol <- rowData(sce)$Symbol
  return(pct_expr_cells_df)
}

# returns dataframe with genes in rows, cell types in cols, and percentage of 
# cells per cell type that express a certain gene for gene cutoff 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# this function was adjusted from Adrien's function interactionmatrix()

#TODO: APPROVED!

extract_matrix <- function(counts, interactions, datasheet, expr_genes){
  # prepared counts matrix from SCE object 
  # interaction database (CellTalkDB)
  # datasheet from prepare_extraction()
  # expr_genes = genes remaining after cutoff (in >X% of cells/identity)
  
  # prepare idents vector to iterate through identities separately by assignment
  cts_emitter <- unfactor(unique(datasheet[datasheet$assignment == "emitter", 
                                           "identity"]))
  cts_receiver <- unfactor(unique(datasheet[datasheet$assignment == "receiver", 
                                            "identity"]))
  
  # empty dataframe of correct dimensions
  interactionmatPlus <- matrix(ncol = ncol(counts), nrow = nrow(interactions))
  interactionmatPlus <- data.frame(interactionmatPlus)
  
  # iterate through all potential interactions i and celltypes j
  for(i in 1:nrow(interactions)){
    for(j in cts_emitter){
      
      #print(interactions$lr_pair[i])
      #print(j)
      
      # if ligand of interaction i is found in the expressed genes of ct j
      if(interactions$ligand_ensembl_gene_id[i] %in% c(expr_genes[[j]])){
        # add to the dataframe
        interactionmatPlus[
          # only into row corresponding to interaction i
          # only into columns that are denoted as cells from identity j
          i, which(datasheet$identity == j)] <- counts[
            # only from the one row in counts
            # whose gene ID is in the receptor col of lrdb_part interaction i
            which(rownames(counts) == interactions[i, "ligand_ensembl_gene_id"]), 
            which(datasheet$identity == j)]
        # fill other slots with NAs
      }else{
        interactionmatPlus[i, which(datasheet$identity == j)] <- NA
      }
    }
    # repeat for receivers
    for(j in cts_receiver){
      if(interactions$ligand_ensembl_gene_id[i] %in% c(expr_genes[[j]])){
        
        interactionmatPlus[i, which(datasheet$identity == j)] <- counts[
          rownames(counts) == interactions[i, "ligand_ensembl_gene_id"], 
          which(datasheet$identity == j)]
        
      }else{
        interactionmatPlus[i, which(datasheet$identity == j)] <- NA
      }
    }
    
    if(i == 1000){print("................ 1000 done")} #track progress 
    
  }
  # remove potential negative values
  interactionmatPlus[interactionmatPlus < 0] <- NA
  return(interactionmatPlus)
}


# returns an interaction matrix with cols = cells and rows = interactions
# extracts counts for ligand genes for emitters and receptor genes for receivers
# extracts counts only for expressed genes that remain after cutoff
















#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# this function was adapted from Adrien's interactionranking() function

interaction_ranking <- function(interactionmatrix_list, 
                                cutoff = TRUE,
                                level = TRUE, 
                                top_level){
  # interactionmatrix_list from interaction_matrix
  # cutoff as indicator whether cutoff was implemented
  # level to indicate if ranks should be leveled to a top_level 
  # top_level = a value for leveling that must be higher than the highest rank

  # Preparation ################################################################
  
  # extract data from list
  interactionmatrix <- interactionmatrix_list[[1]]
  datasheet <- interactionmatrix_list[[2]]
  
  # reformat datasheet into two
  datasheetemi <- datasheet[datasheet[,3] %in% "emitter",]
  datasheetrec <- datasheet[datasheet[,3] %in% "receiver",]
  
  # empty objects 
  name <- vector()
  mat <- matrix(
    ncol = length(unique(datasheetemi[,2])) * length(unique(datasheetrec[,2])),
    nrow = nrow(interactionmatrix))
  Score <- mat
  rankLigs <- mat
  rankRecs <- mat
  
  # 1 will be added during each iteration 
  l <- 1
  
  # rank the ligand and receptor genes #########################################
  
  # for each ligand i and for each receptor j
  for(i in 1:length(unique(datasheetemi[,2]))){
    for (j in 1:length(unique(datasheetrec[,2]))){
      # get all positions of the current ligand i 
      pos_lig <- which(datasheet[,2] %in% unique(datasheetemi[,2])[i])
      # rank the mean expression values for those positions
      if(length(pos_lig) > 1){
        rankLig <- dense_rank(rowMeans(
          interactionmatrix[,pos_lig]))
      }else if(length(pos_lig) == 1){
        # if there is only one cell = one column, rowMeans does not work
        rankLig <- dense_rank(interactionmatrix[,pos_lig])
      }
      
      # same for receptors
      pos_rec <- which(datasheet[,2] %in% unique(datasheetrec[,2])[j])
      if(length(pos_rec) > 1){
        rankRec <- dense_rank(rowMeans(
          interactionmatrix[,pos_rec]))
      }else if(length(pos_rec) == 1){
        rankRec <- dense_rank(interactionmatrix[,pos_rec])
      }
      
      # level the ranks ########################################################
      
      # level ranks by adding up to top_level to reduce influences on score
      if(level == TRUE && cutoff == FALSE){
        # a rank of 1 indicates 0 expression (no cut off)

        # determine the highest rank
        max_ints_lig <- max(rankLig, na.rm = TRUE)
        max_ints_rec <- max(rankRec, na.rm = TRUE)
        # determine the resulting lowest rank
        min_rank_lig <- top_level-max_ints_lig
        min_rank_rec <- top_level-max_ints_rec
        # stop the function if the resulting lowest rank is <0
        if(min_rank_lig < 0 | min_rank_rec < 0){
          stop("top_level is too low, set top_level higher")
        }
        
        # add resulting lowest rank to remaining ranks to match the highest rank 
        # but only if ranks are not 1, which indicates expression of 0
        rankLig[rankLig != 1] <- rankLig[rankLig != 1]+min_rank_lig
        rankRec[rankRec != 1] <- rankRec[rankRec != 1]+min_rank_rec
        # remove 1s to avoid tiny scores after score normalization later
        rankLig[rankLig == 1] <- 0
        rankRec[rankRec == 1] <- 0      
        
      }else if(level == TRUE && cutoff == TRUE){
        # a rank of 1 indicates lowest expressed rank (no non-expressing cells)

        max_ints_lig <- max(rankLig, na.rm = TRUE)
        max_ints_rec <- max(rankRec, na.rm = TRUE)
        min_rank_lig <- top_level-max_ints_lig
        min_rank_rec <- top_level-max_ints_rec
        
        if(min_rank_lig < 0 | min_rank_rec < 0){
          stop("max rank level is too low, set max rank level higher")
        }
        
        # add lowest rank only if ranks are not NA, indicating expression of 0
        rankLig[!is.na(rankLig)] <- rankLig[!is.na(rankLig)]+min_rank_lig
        rankRec[!is.na(rankRec)] <- rankRec[!is.na(rankRec)]+min_rank_rec
      }
    
      # put in empty mat, with col corresponding to the number of iteration
      rankLigs[,l] <- rankLig
      rankRecs[,l] <- rankRec

      
      # calculate Rank Sums, Score and match cell types ########################
      
      Score[,l] <- rankLig + rankRec
      
      # store the colnames for Score for later
      name[l] <- paste0(unique(datasheetemi[,2])[i], "&",
                        unique(datasheetrec[,2])[j])
      l=l+1
    }
  }
  
  # Score normalization for each pair of cells 
  for (i in 1:ncol(Score)){
    Score[,i] <- Score[,i]/max(Score[,i], na.rm = TRUE)
  }
  colnames(Score) <- name
  
  interactionranking_list <- list()
  interactionranking_list[[1]] <- Score
  interactionranking_list[[2]] <- rankRecs
  interactionranking_list[[3]] <- rankLigs
  interactionranking_list[[4]] <- datasheet  # keep datasheet
  return(interactionranking_list)
}

# interaction_ranking
# returns a list of three matrices containing for all ctp and interactions
# - interaction score 
# - receptor ranks 
# - ligand ranks 
# - and the datasheet


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


interaction_list <- function(interactionranking_list, lrdb_comp){
  # interactionranking_list from interaction_ranking
  # lrdb_comp containing more info on all interactions
  
  # extract info 
  ir_matrix <- list(interactionranking_list[[1]], # Score
                    interactionranking_list[[2]], # Receptorranks
                    interactionranking_list[[3]]) # Ligandranks
  datasheet <- interactionranking_list[[4]]
  lrdb <- lrdb_comp
  
  # add names 
  names(ir_matrix) <- c("Score", "Receptorrank", "Ligandrank")
  rownames(ir_matrix$Score) <- lrdb$interactions$interaction_pair
  rownames(ir_matrix$Receptorrank) <- lrdb$interactions$interaction_pair 
  rownames(ir_matrix$Ligandrank) <-lrdb$interactions$interaction_pair
  colnames(ir_matrix$Receptorrank) <- colnames(ir_matrix$Score)
  colnames(ir_matrix$Ligandrank) <- colnames(ir_matrix$Score)
  
  # coerce into dataframes
  ir_matrix$Score <- as.data.frame(ir_matrix$Score)
  ir_matrix$Receptorrank <- as.data.frame(ir_matrix$Receptorrank)
  ir_matrix$Ligandrank <- as.data.frame(ir_matrix$Ligandrank)
  
  # add colData-like slot = Celltypes ##########################################
  ir_matrix$Celltypes <- data.frame(
    ct_pair = colnames(ir_matrix$Score),
    emitter = vector(mode = "character", length = ncol(ir_matrix$Score)),
    receiver = vector(mode = "character", length = ncol(ir_matrix$Score))
  )
  
  # add annotation as emitter or receiver for better overview
  ct <- unique(datasheet$celltype)
  
  for(i in ct){
    if(datasheet$annotation[datasheet$celltype == i][1] == "emitter"){
      ir_matrix$Celltypes$emitter[grep(i, ir_matrix$Celltypes$ct_pair)] <- i
    }else if(datasheet$annotation[datasheet$celltype == i][1] == "receiver"){
      ir_matrix$Celltypes$receiver[grep(i, ir_matrix$Celltypes$ct_pair)] <- i
    }
  }
  
  # add rowData-like slot = Interactions #######################################
  ir_matrix$Interactions <- data.frame(
    interaction_pair = lrdb$interactions$interaction_pair,
    ligand_symbol = lrdb$interactions$ligand_symbol,
    ligand_ensembl_id = lrdb$interactions$ligand_ensembl_id,
    receptor_symbol = lrdb$interactions$receptor_symbol,
    receptor_ensembl_id = lrdb$interactions$receptor_ensembl_id,
    interaction_type = lrdb$interactions$interaction_type
  )
  return(ir_matrix)
}

# transforms a list from interaction_ranking() into into another list
# with more information 
