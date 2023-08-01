# functions for calculation of cell type interactomes 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

## PERCENTAGE OF EXPRESSED CELLS (expr_cells_perc)

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

# EXTRACT MATRIX
# this function was adjusted from Adrien's function interactionmatrix()

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
    
    rownames(interactionmatPlus)[i] <- interactions$lr_pair[i]
    print(interactions$lr_pair[i])
    
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
            # whose gene ID is in the ligand col of lrdb_part interaction i
            which(rownames(counts) == interactions[i, "ligand_ensembl_gene_id"]), 
            which(datasheet$identity == j)]
        # fill other slots with NAs
      }else{
        interactionmatPlus[i, which(datasheet$identity == j)] <- NA
      }
    }
    # repeat for receivers
    for(j in cts_receiver){
      # if receptor of interaction i is found in the expressed genes of ct j
      if(interactions$receptor_ensembl_gene_id[i] %in% c(expr_genes[[j]])){
        
        # add to the dataframe only into row corresponding to interaction i
        # only into columns that are denoted as cells from identity j
        interactionmatPlus[i, which(datasheet$identity == j)] <- counts[
          # only from the one row in counts
          # whose gene ID is in the receptor col of lrdb_part interaction i
          rownames(counts) == interactions[i, "receptor_ensembl_gene_id"], 
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

# INTERACTION RANKING
# this function was adapted from Adrien's interactionranking() function

# TODO: Shorten once testing has been established
interaction_ranking <- function(interaction_mat,
                                datasheet,
                                cutoff = TRUE,
                                level = TRUE, 
                                top_level){
  # interactionmatrix_list from interaction_matrix
  # cutoff as indicator whether cutoff was implemented
  # level to indicate if ranks should be leveled to a top_level 
  # top_level = a value for leveling that must be higher than the highest rank
  
  # Preparation ################################################################
  
  # reformat datasheet into two
  datasheetemi <- datasheet[datasheet$assignment == "emitter",]
  datasheetrec <- datasheet[datasheet$assignment == "receiver",]
  
  # generate empty objects 
  name <- vector()
  mat <- matrix(
    ncol = length(unique(datasheetemi$identity)) * 
      length(unique(datasheetrec$identity)),
    nrow = nrow(interaction_mat))
  Score <- mat
  rankLigs <- mat
  rankRecs <- mat
  
  # 1 will be added during each iteration 
  l <- 1
  
  # rank the ligand and receptor genes #########################################
  
  # for each emitter i and for each receiver j
  for(i in 1:length(unique(datasheetemi$identity))){
    for (j in 1:length(unique(datasheetrec$identity))){
      # get all positions of the current receiver cells i in the matrix
      pos_lig <- which(datasheet$identity %in% unique(datasheetemi$identity)[i])
      # rank the mean expression values of all genes for those positions
      # due to prior construction, there are only ligand genes in emitters slots
      datasheet[pos_lig,]
      interaction_mat[,pos_lig]
      if(length(pos_lig) > 1){
        rankLig <- dense_rank(rowMeans(
          interaction_mat[,pos_lig]))
      }else if(length(pos_lig) == 1){
        # if there is only one cell = one column, rowMeans does not work
        # this should only be relevant during testing
        rankLig <- dense_rank(interaction_mat[,pos_lig])
      }
      table(rankLig)
      
      # same for receptors
      pos_rec <- which(datasheet$identity %in% unique(datasheetrec$identity)[j])
      datasheet[pos_rec,]
      interaction_mat[,pos_rec]
      
      if(length(pos_rec) > 1){
        rankRec <- dense_rank(rowMeans(
          interaction_mat[,pos_rec]))
      }else if(length(pos_rec) == 1){
        rankRec <- dense_rank(interaction_mat[,pos_rec])
      }
      table(rankRec)
      
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
      name[l] <- paste0(unique(datasheetemi$identity)[i], "&",
                        unique(datasheetrec$identity)[j])
      l <- l+1
    }
  }
  
  # Score normalization for each pair of cells 
  for (i in 1:ncol(Score)){
    Score[,i] <- Score[,i]/max(Score[,i], na.rm = TRUE)
  }
  print(name)
  colnames(Score) <- name
  
  interactionranking_list <- list()
  interactionranking_list[["Score"]] <- Score
  interactionranking_list[["Ligandrank"]] <- rankLigs
  interactionranking_list[["Receptorrank"]] <- rankRecs
  interactionranking_list[["datasheet"]] <- datasheet  # keep datasheet
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

# MAKE AN INTERACTION LIST

make_interaction_list <- function(interaction_ranking, lrdb){
  # interactionranking_list from interaction_ranking
  # lrdb_comp containing more info on all interactions
  
  # extract info 
  interaction_list <- list(interaction_ranking$Score, 
                           interaction_ranking$Receptorrank,
                           interaction_ranking$Ligandrank) 
  datasheet <- interaction_ranking$datasheet
  
  # add names 
  # the sequence of rows is the same is interaction_mat is the same as lrdb
  names(interaction_list) <- c("Score", "Receptorrank", "Ligandrank")
  rownames(interaction_list$Score) <- lrdb$lr_pair
  rownames(interaction_list$Receptorrank) <- lrdb$lr_pair
  rownames(interaction_list$Ligandrank) <- lrdb$lr_pair
  colnames(interaction_list$Receptorrank) <- colnames(interaction_list$Score)
  colnames(interaction_list$Ligandrank) <- colnames(interaction_list$Score)
  
  # coerce into dataframes
  interaction_list$Score <- as.data.frame(interaction_list$Score)
  interaction_list$Receptorrank <- as.data.frame(interaction_list$Receptorrank)
  interaction_list$Ligandrank <- as.data.frame(interaction_list$Ligandrank)
  
  # add colData-like slot = Celltypes ##########################################
  interaction_list$Identities <- data.frame(
    ident_pair = colnames(interaction_list$Score),
    emitter = vector(mode = "character", length = ncol(interaction_list$Score)),
    receiver = vector(mode = "character", length = ncol(interaction_list$Score))
  )
  
  # add annotation as emitter or receiver for better overview
  idents <- unique(datasheet$identity)
  
  for(i in idents){
    if(datasheet$assignment[datasheet$identity == i][1] == "emitter"){
      interaction_list$Identities$emitter[
        grep(i, interaction_list$Identities$ident_pair)] <- i
    }else if(datasheet$assignment[datasheet$identity == i][1] == "receiver"){
      interaction_list$Identities$receiver[
        grep(i, interaction_list$Identities$ident_pair)] <- i
    }
  }
  
  # add rowData-like slot = Interactions #######################################
  interaction_list$Interactions <- data.frame(
    interaction_pair = lrdb$lr_pair,
    ligand_symbol = lrdb$ligand_gene_symbol,
    ligand_ensembl_id = lrdb$ligand_ensembl_gene_id,
    receptor_symbol = lrdb$receptor_gene_symbol,
    receptor_ensembl_id = lrdb$receptor_ensembl_gene_id
  )
  return(interaction_list)
}

# transforms a list from interaction_ranking() into into another list
# with more information 

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# EXTRACT INFO PER IDENTITY PAIR
extract_ident_pair_info <- function(int_list){
  #  interactome object
  
  # empty objects
  return_list <- list()
  return_df <- data.frame(
    "ident_pairs" = colnames(int_list$Score)
  )
  return_df$nr_ints <- vector(length = nrow(return_df))
  return_df$emitter <- vector(length = nrow(return_df))
  return_df$receiver <- vector(length = nrow(return_df))
  
  # for each ctp in the cpi object, perform separately
  # this will yield ncol dataframes per cpi object in a large return_list
  for(i in 1:ncol(int_list$Score)){
    temp_df <- data.frame(
      "ident_pair" = rep(colnames(int_list$Score)[i], nrow(int_list$Score)),
      "interaction" = rownames(int_list$Score),
      "score" = int_list$Score[,i],
      "emitter" = rep(int_list$Identities$emitter[i], nrow(int_list$Score)),
      "ligandrank" = int_list$Ligandrank[,i],
      "receiver" = rep(int_list$Identities$receiver[i], nrow(int_list$Score)),
      "receptorrank" = int_list$Receptorrank[,i]
    )
    
    # remove all non-detected interactions, then order from high to low score
    temp_df <- temp_df[!is.na(temp_df$score),]
    temp_df <- temp_df[order(temp_df$score, decreasing = TRUE),]
    return_list[[i]] <- temp_df
    
    # enter one row per ctp into a df collecting the nr of interactions per ctp
    # for overview
    return_df$nr_ints[i] <- nrow(temp_df) 
    return_df$emitter[i] <- temp_df$emitter[1]
    return_df$receiver[i] <- temp_df$receiver[1]
  }
  names(return_list) <- colnames(int_list$Score)
  return_df <- return_df[order(return_df$nr_ints, decreasing = TRUE),]
  print(return_df)
  
  # the resulting object contains an overview df of all ctp interactions, 
  # then separate dfs for each separate ctp
  return(list("overview" = return_df, 
              "separate" = return_list))
}

# this function extracts important info on interactions at IDENTITIY PAIR level
# can only be used if cutoff = TRUE and level = TRUE and no subsetting required
# otherwise, use extract_ident_pair_info_raw()

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# EXTRACT INFO PER IDENTITY

extract_ident_info <- function(ipi_list){
  
  # put all info from ctpints lists into one big dataframe
  ipi_total <- ipi_list[[2]][[1]]
  for(i in 2:length(ipi_list[[2]])){
    ipi_total <- rbind(ipi_total, ipi_list[[2]][[i]])
  }
  
  # then rearrange by cell types, first emitters
  nr_emitters <- length(unique(ipi_total$emitter))
  nr_receivers <- length(unique(ipi_total$receiver))
  
  ipi_total_list <- list()
  for(i in 1:nr_emitters){
    ipi_total_list[[i]] <- ipi_total[
      ipi_total$emitter == unique(ipi_total$emitter)[i],]
    ipi_total_list[[i]]$identity <- unique(ipi_total$emitter)[i]
    ipi_total_list[[i]]$assignment <- "emitter"
  }
  
  # then receivers
  for(i in 1:nr_receivers){
    
    j <- nr_emitters + i
    ipi_total_list[[j]] <- ipi_total[
      ipi_total$receiver == unique(ipi_total$receiver)[i],]
    ipi_total_list[[j]]$identity <- unique(ipi_total$receiver)[i]
    ipi_total_list[[j]]$assignment <- "receiver"
  }
  # ipi_total_list now contains ident_pair_infos for each emitter and receiver
  
  # now use this data to make a ctints_list with all relevant information
  ident_info_list <- lapply(ipi_total_list, function(ipi_total){
    
    ident_info <- data.frame(
      interactions = unique(ipi_total$interaction),
      identity = rep(ipi_total$identity, 
                     length = length(unique(ipi_total$interaction))),
      assignment = rep(ipi_total$assignment,
                       length = length(unique(ipi_total$interaction)))
    )
    
    # add the relevant rank
    if(ident_info$assignment[1] == "emitter"){
      ident_info$rank <-  ipi_total$ligandrank[
        match(ident_info$interactions, ipi_total$interaction)]
    }else if(ident_info$assignment[1] == "receiver"){
      ident_info$rank <-  ipi_total$receptorrank[
        match(ident_info$interactions, ipi_total$interaction)]
    }
    
    ident_info <- ident_info[order(ident_info$rank, decreasing = TRUE),]
    return(ident_info)
  })
  
  names(ident_info_list) <- c(unique(ipi_total$emitter), 
                              unique(ipi_total$receiver))
  
  # make an overview df
  overview_df <- data.frame(
    identity =  names(ident_info_list)
  )
  
  for(i in 1:length(ident_info_list)){
    overview_df$nr_ints[i] <- nrow(ident_info_list[[i]])
    overview_df$assignment[i] <- ident_info_list[[i]]$assignment[1]
  }
  overview_df <- overview_df[order(overview_df$nr_ints, decreasing = TRUE),]
  print(overview_df)
  
  return(list(overview = overview_df, 
              separate = ident_info_list))
  
}

# this function uses the ctpints list to extract info on identity level

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# EXTRACT LIGAND/RECEPTOR INFO PER IDENTITY

extract_lrs_info <- function(idi_list){
  # object with info on interactions on CELL TYPE level from extract_ct_info()
  
  if(idi_list$assignment[1] == "emitter"){
    # remove the receiver from the lri name to obtain only the ligand
    idi_list$interactions <- gsub("[_][[:alnum:]]+", "", idi_list$interactions)
    ligands <- unique(idi_list$interactions)
    
    return_df <- data.frame(
      identity = rep(idi_list$identity[1], length(ligands)),
      assignment = rep("emitter", length(ligands)),
      sigmols = ligands
    )
  }
  if(idi_list$assignment[1] == "receiver"){
    # remove the emitter from the lri name to obtain only the receptor
    idi_list$interactions <- gsub("[[:alnum:]]+[_]", "", idi_list$interactions)
    receptors <- unique(idi_list$interactions)
    
    return_df <- data.frame(
      identity = rep(idi_list$identity[1], length(receptors)),
      assignment = rep("receiver", length(receptors)),
      sigmols = receptors
    )
  }
  #print(return_df)
  return(return_df)
}

# for each ct = element in ctints_list[[2]], the expressed l or r are obtained
# returns a list with identical architecture to ctints_list[[2]] 

#-------------------------------------------------------------------------------

extract_lrs_nrs <- function(int_list, lrs_list){
  # int_list = cell type pair interactome object
  # lrs_list = info ligands or receptors per identity from extract_lrs_info()
  
  print_list <- list()
  
  emitter_ct <- unique(int_list$Identities$emitter)
  receiver_ct <- unique(int_list$Identities$receiver)
  idents <- c(emitter_ct, receiver_ct)
  
  print_df <- data.frame(
    identity = vector(mode = "character", length = length(idents)),
    nr_lrs = vector(mode = "numeric", length = length(idents)),
    assignment = vector(mode = "character", length = length(idents))
  )
  
  for(i in 1:length(idents)){
    print(i)
    print(idents[i])
    print(lrs_list[[i]])
    print_df$identity[i] <- lrs_list[[i]]$identity[i]
    print_df$nr_lrs[i] <- nrow(lrs_list[[i]])
    print_df$assignment[i] <- lrs_list[[i]]$assignment[i]
  }
  
  print_df <- print_df[order(print_df$nr_lrs, decreasing = TRUE),]
  
  return(print_df)
}

# returns a df per cpi object with the nr of l or r per cell type
# corresponds to "overview" df of extract_ctp_info() and extract_ct_info()
