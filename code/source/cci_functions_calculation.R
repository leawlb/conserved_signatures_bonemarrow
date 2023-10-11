
library(tidyverse)
library(dplyr)
library(SingleCellExperiment)

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# functions for calculation of cell type interactomes == CCI objects


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## MAIN FUNCTION: STEP 1 
## extract the percentage of cells per identity that express a given gene:
## GET PERCENTAGE OF EXPRESSING CELLS = perc_expr_cells

# returns dataframe with gene IDs = rows, cell types = cols
# with percentage of cells per identity that express a given gene 

perc_expr_cells <- function(sce, assay_use = NULL){
  
  # purpose: extract % of cells per identity that express a given gene
  
  # pass through parameters
  sce <- sce                # SCE object, should be pre-processed
  assay_use <- assay_use    # assay to use, "downsampled" by default
  if(is.null(assay_use) == TRUE){assay_use <- "downsampled"}
  
  # for each identity i
  ident_list <- as.list(levels(sce$Identity))
  perc_list <- lapply(ident_list, function(i){
    
    # subset sce by identity and into assay to use
    sce_temp <- sce[,colData(sce)$Identity == i] 
    assay <- assays(sce_temp)[[which(names(assays(sce_temp)) == assay_use)]]
    #track progress
    print(paste("........ ", i, ": ", ncol(sce_temp), " cells", sep = ""))
    
    # for each gene g
    gene_list <- as.list(rowData(sce)$Symbol)
    perc_exp_list <- lapply(gene_list, function(g){
      
      # % of cells per ct that have at least 1 counts for gene g in identity i
      nr_cell_expr <- length(which(assay[g,] > 0))
      perc_exp <- nr_cell_expr/ncol(sce_temp)
      return(perc_exp)
    })
    
    perc_exp <- unlist(perc_exp_list)
    return(perc_exp)
  })
  
  # convert into one df
  names(perc_list) <- levels(sce$Identity)
  perc_df <- as.data.frame(bind_cols(perc_list))
  # add data on genes
  rownames(perc_df) <- rowData(sce)$ENSMUS_ID
  perc_df$gene_symbol <- rowData(sce)$Symbol
  
  return(perc_df)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## HELPER FUNCTION 
## remove genes that are not expressed in at least "cutoff" percent of cells:
## CUTOFF OF NON-EXPRESSED GENES = nexpr_gene_cutoff

# returns list of genes per identity that are expressed above "cutoff" percent
# of cells per identity, for subsetting

nexpr_gene_cutoff <- function(ident, perc_df, cutoff){
  
  # purpose: keep genes that are expressed in at least "cutoff" % of cells:
  
  ident <- ident       # current identity
  perc_df <- perc_df   # df with % of expressing cells from perc_expr_cells()
  cutoff <- cutoff     # cutoff value between 0 and 1 (e.g. 0.2 = 20%)

  # only keep genes above cutoff
  expr_genes <- rownames(perc_df)[
    which(perc_df[,colnames(perc_df) == ident] >= cutoff)]
  
  return(expr_genes)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## HELPER FUNCTION 
## various small steps to prepare for extraction of the interaction matrix:
## PREPARE EXTRACTION

# returns list of various objects required for downstream calculation

prepare_extraction <- function(sce, assay_use){
  
  # purpose: prepare a list of objects required for downstream calculation
  
  sce <- sce               # SCE object, should be pre-processed
  assay_use <- assay_use   # assay to use, "downsampled" by default
  if(is.null(assay_use) == TRUE){assay_use <- "downsampled"}
  
  # create countsmatrix and datasheet for extract_matrix() 
  countsmatrix <- assays(sce)[[which(names(assays(sce)) == assay_use)]]
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





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## MAIN FUNCTION: STEP 2 EXTRACT_INTERACTION_MATRIX
## extract relevant gene expression values into interaction matrix
## ligand expression for emitters, receptor expression for receivers
## EXTRACT INTERACTION MATRIX = extract_matrix()

# this function was adjusted from Adrien Jolly's function interactionmatrix()
# https://github.com/AdrienJolly/CellInteractionScores/blob/master/InteractionRanking.r

# returns an interaction matrix as dataframe
# cols = cells and rows = ligand receptor interactions
# with ligand gene counts for emitters and receptor genes counts for receivers
# only for relevant expressed genes 

extract_matrix <- function(counts, datasheet, lrdb, expr_genes){
  
  # purpose: extract relevant gene expression values into interaction matrix

  # pass through parameters
  counts <- counts          # prepared counts matrix, from prepare_extraction()
  datasheet <- datasheet    # df with info on cells, from prepare_extraction()
  lrdb <- lrdb              # ligand receptor database 
  expr_genes <- expr_genes  # genes after cutoff from nexpr_gene_cutoff()
  
  # construct empty df to iterate through to fill
  int_df <- data.frame(matrix(ncol = ncol(counts), nrow = nrow(lrdb)))
  colnames(int_df) <- colnames(counts)
  rownames(int_df) <- lrdb$lr_pair
  
  # ct vectors to iterate through
  cts_emitter <- unfactor(unique(datasheet[datasheet$assignment == "emitter", 
                                           "identity"]))
  cts_receiver <- unfactor(unique(datasheet[datasheet$assignment == "receiver", 
                                            "identity"]))
  
  print(cts_emitter)
  print(cts_receiver)
  
  # transpose int_df for use with mapply
  mint_df <- as.data.frame(t(int_df))
  # carry colname = interaction i into mapply 
  mint_df[1,] <- colnames(mint_df)
  
  # add expr values from counts into empty interaction matrix (int_df)
  add_counts <- function(x){
    
    x <- x
    # extract interaction i, turn slot into NA and num again
    int_curr <- x[1]
    x[1] <- NA
    x <- as.numeric(x)

    # get the ligand and receptor genes of current interaction i
    lig_curr <- lrdb$ligand_ensembl_gene_id[lrdb$lr_pair == int_curr]
    rec_curr <- lrdb$receptor_ensembl_gene_id[lrdb$lr_pair == int_curr]
    
    # get the ligand genes expressed in emitter cell types
    for(e in cts_emitter){
      
      # only if ligand gene is in "expressed genes" that made the cut-off
      if(lig_curr %in% c(expr_genes[[e]])){
        
        # positions of the cells in counts cols, and int_df rows = pos in x
        emi_pos <- which(datasheet$identity == e)
        
        # extract from counts: row corresponding to lig_gene_curr,
        # col corresponding to cells from emitter e, and
        # transfer into emi_df: 
        # emi_df only has 1 row corresponding to interaction nr i,
        # only has cols corresponding to the same cells from emitter e
        x[emi_pos] <- counts[which(rownames(counts) == lig_curr), emi_pos]
      }
    }
    
    # same for receivers
    for(r in cts_receiver){
      if(rec_curr %in% c(expr_genes[[r]])){
        rec_pos <- which(datasheet$identity == r)
        x[rec_pos] <- counts[which(rownames(counts) == rec_curr), rec_pos]
      }
    }
    return(x)
  }
  
  # transpose back, adjust to required output as dataframe
  mint_df_out <- mapply(x = mint_df, FUN = add_counts)
  int_df_return <-t(mint_df_out)
  colnames(int_df_return) <- colnames(int_df)
  
  stopifnot(length(which(!is.na(int_df_return))) > 0)
  stopifnot(!unname(rowSums(int_df_return, na.rm = TRUE)) < 0)
  stopifnot(identical(colnames(int_df_return), datasheet$sample))
  stopifnot(identical(rownames(int_df_return), lrdb$lr_pair))
  
  int_df_return <- as.data.frame(int_df_return)
  return(int_df_return)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## MAIN FUNCTION: STEP 3 
## rank ligand and receptor expression, score resulting interactions:
## INTERACTION SCORING

# this function was adapted from Adrien Jolly's interactionranking() function
# https://github.com/AdrienJolly/CellInteractionScores/blob/master/InteractionRanking.r

# returns a list of three matrices containing for all ctp and interactions
# - interaction score (ints = rows, identiy pairs = cols)
# - receptor ranks 
# - ligand ranks 
# - and the datasheet

# TODO: Shorten once testing has been established
interaction_scoring_testing <- function(interaction_mat,
                                        datasheet,
                                        cutoff = TRUE,
                                        level = TRUE, 
                                        top_level){
  
  # interaction_mat = df of interaction matrix obtained from extract_matrix()
  # datasheet = datasheet with info on cells, from prepare_extraction()
  # cutoff = set true if cutoff is implemented
  # level = set true if top_level is implemented
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
  
  # Score normalization for each pair of identities 
  for (i in 1:ncol(Score)){
    # set scores so that the lowest Score is 1
    Score[,i] <- Score[,i] - (min(Score[,i], na.rm = TRUE) -1) # remove if it doesn't work
    Score[,i] <- Score[,i]/max(Score[,i], na.rm = TRUE)
  }
  
  # instead of normalizing to max, which results in loss of values between
  # 0 and ~0.6, use percent_rank to obtain the percentile of each score within
  # all scores
  #for (i in 1:ncol(Score)){
  #  Score[,i] <- percent_rank(Score[,i])
  #}
  #print(name)
  colnames(Score) <- name
  rownames(Score) <- rownames(interaction_mat)
  
  colnames(rankLigs) <- name
  rownames(rankLigs) <- rownames(Score)
  
  colnames(rankRecs) <- name
  rownames(rankRecs) <- rownames(Score)
  
  
  interactionscoring_list <- list()
  interactionscoring_list[["Score"]] <- Score
  interactionscoring_list[["Ligandrank"]] <- rankLigs
  interactionscoring_list[["Receptorrank"]] <- rankRecs
  interactionscoring_list[["datasheet"]] <- datasheet  # keep datasheet
  return(interactionscoring_list)
}

#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------


## MAIN FUNCTION: STEP 3 INTERACTION_SCORING
## rank ligand and receptor expression, then score resulting interactions:
## INTERACTION SCORING

# this function was adapted from Adrien's interactionranking() function

# returns a list of three matrices containing for all ctp and interactions
# - interaction score 
# - receptor ranks 
# - ligand ranks 
# - and the datasheet

# just slightly quicker
interaction_scoring <- function(int_df,
                                datasheet,
                                top_level){
  
  # purpose: rank ligand and receptor expression, score interactions
  
  int_df <- int_df        # df of interaction matrix from extract_matrix()
  datasheet <- datasheet  # df with info on cells, from prepare_extraction()
  top_level <- top_level  # for leveling, must be higher than the highest rank

  # Preparation ################################################################
  
  # reformat datasheet, get cell type lists per assignment
  datasheetemi <- datasheet[datasheet$assignment == "emitter",]
  datasheetrec <- datasheet[datasheet$assignment == "receiver",]
  emi_list <- as.list(unfactor(unique(datasheetemi$identity)))
  rec_list <- as.list(unfactor(unique(datasheetrec$identity)))
  
  # generate empty mat of correct structure for mapply 
  mat <- data.frame(
    matrix(ncol = length(unfactor(unique(datasheetemi$identity))) * 
             length(unfactor(unique(datasheetrec$identity))), 
           nrow = nrow(int_df)))
  
  # construct names of identity pairs
  names_list <- lapply(emi_list, function(e){
    names <- lapply(rec_list, function(r){
      name <- paste0(e, "&", r)
      return(name)
    })
    names <- unlist(names)
    return(names)
  })
  names_mat <- unlist(names_list)
  colnames(mat) <- names_mat
  rownames(mat) <- rownames(int_df)
  
  # carry info on identity pairs into mapply
  mat[1,] <- colnames(mat)
  
  # rank ligands
  rank_ligs <- function(x){
    
    # extract interaction i, turn slot into NA again
    idp_curr <- x[1]
    x[1] <- NA
    x <- as.numeric(x)
    
    # get positions of identity e in mat (same as datasheet)
    emi <- str_split(idp_curr, "&")[[1]][1]
    pos_emi <- which(datasheet$identity == emi)
    
    # rank all genes based on their mean expression values for those positions
    # due to prior df construction there are only ligand genes in emi slots
    rankLig <- dense_rank(rowMeans(int_df[,pos_emi]))
    
    # level ranks 
    # determine the highest rank
    max_ints_lig <- max(rankLig, na.rm = TRUE)
    # determine the resulting lowest rank
    min_rank_lig <- top_level-max_ints_lig
          
    if(min_rank_lig < 0){stop("min rank level < 0")}
    rankLig[!is.na(rankLig)] <- rankLig[!is.na(rankLig)]+min_rank_lig
    
    
    return(rankLig)
  }
  
  # rank receptors and level to top_level for better comparability
  rank_recs <- function(x){
    
    idp_curr <- x[1]
    x[1] <- NA
    x <- as.numeric(x)
    
    rec <- str_split(idp_curr, "&")[[1]][2]
    pos_rec <- which(datasheet$identity == rec)
    
    rankRec <- dense_rank(rowMeans(int_df[,pos_rec]))
    
    # level ranks
    max_ints_rec <- max(rankRec, na.rm = TRUE)
    min_rank_rec <- top_level-max_ints_rec
    if(min_rank_rec < 0){stop("min rank level < 0")}
    rankRec[!is.na(rankRec)] <- rankRec[!is.na(rankRec)]+min_rank_rec
    
    return(rankRec)
  }
  rankLigs <- mapply(mat, FUN = rank_ligs)
  rankRecs <- mapply(mat, FUN = rank_recs)
  
  rownames(rankLigs) <- rownames(int_df)
  rownames(rankRecs) <- rownames(int_df)
  colnames(rankLigs) <- names_mat
  colnames(rankRecs) <- names_mat
  
  calc_score <- function(x = mat){
    
    idp_curr <- x[1]
    x[1] <- NA
    x <- as.numeric(x)
    
    rankLigs_curr <- rankLigs[,colnames(rankLigs) == idp_curr]
    rankRecs_curr <- rankRecs[,colnames(rankRecs) == idp_curr]
    
    score <-  rankLigs_curr + rankRecs_curr
    
    # adjust Scores so that lowest score is 1
    score <- score - (min(score, na.rm = TRUE) -1) 
    # normalize, resulting in scores from almost 0 to 1
    score <- score/max(score, na.rm = TRUE)
    
    return(score)
  }
  
  Score <- mapply(mat, FUN = calc_score)
  rownames(Score) <- rownames(int_df)
  colnames(Score) <- names_mat
  
  int_scoring_list <- list()
  int_scoring_list[["Score"]] <- Score
  int_scoring_list[["Ligandrank"]] <- rankLigs
  int_scoring_list[["Receptorrank"]] <- rankRecs
  int_scoring_list[["datasheet"]] <- datasheet  # keep datasheet
  return(int_scoring_list)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## MAIN FUNCTION: STEP 4 
## re-organize list of scored interactions into more convenient format:
## MAKE INTERACTION LIST

# returns a CCI object = a list of various dataframes on scores, ranks, 
# identities, and interactions

make_interaction_list <- function(interaction_scores, lrdb){
  
  # purpose: make final CCI object, reorganize interaction scores list
  
  interaction_scores <- interaction_scores  # list from interaction_scoring()
  lrdb <- lrdb                              # ligand receptor database 
  
  # extract info 
  interaction_list <- list(interaction_scores$Score, 
                           interaction_scores$Receptorrank,
                           interaction_scores$Ligandrank) 
  datasheet <- interaction_scores$datasheet
  
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

# from here on, interaction_list will be called CCI object 
