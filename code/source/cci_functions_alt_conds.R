

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------


## ALTERNATIVE FUNCTION: STEP 3 INTERACTION_SCORING
## rank ligand and receptor expression, then score resulting interactions:
## INTERACTION SCORING

# this function was adapted from Adrien's interactionranking() function
# this function is slightly changed from interaction_scoring:
# - added option for no_level (no levelling)
# - added option for no_cutoff and no_preprocessing (adding NAs)

# returns a list of three matrices containing for all ctp and interactions
# - interaction score 
# - receptor ranks 
# - ligand ranks 
# - and the datasheet

# just slightly quicker
interaction_scoring_alt <- function(int_df,
                                    datasheet,
                                    top_level){
  
  # purpose: rank ligand and receptor expression, score interactions
  
  int_df <- int_df        # df of interaction matrix from extract_matrix()
  datasheet <- datasheet  # df with info on cells, from prepare_extraction()
  top_level <- top_level  # for leveling, must be higher than the highest rank
  
  # Preparation ################################################################
  
  # reformat datasheet into two, get cell type lists per assignment
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
    
    # rank the mean expression values of all genes for those positions
    # due to prior df construction there are only ligand genes in emi slots
    rankLig <- dense_rank(rowMeans(int_df[,pos_emi]))
    
    # level ranks (except for testing no_level)
    if(!is.null(top_level)){
      # determine the highest rank
      max_ints_lig <- max(rankLig, na.rm = TRUE)
      # determine the resulting lowest rank
      min_rank_lig <- top_level-max_ints_lig
      
      if(min_rank_lig < 0){stop("min rank level < 0")}
      rankLig[!is.na(rankLig)] <- rankLig[!is.na(rankLig)]+min_rank_lig
    }
    
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
    
    # level ranks (except for testing no_level)
    if(!is.null(top_level)){
      max_ints_rec <- max(rankRec, na.rm = TRUE)
      min_rank_rec <- top_level-max_ints_rec
      if(min_rank_rec < 0){stop("min rank level < 0")}
      rankRec[!is.na(rankRec)] <- rankRec[!is.na(rankRec)]+min_rank_rec
    }
    
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
    
    # if cutoff was performed, there are NAs for non-expressed genes
    # if no cutoff was performed, there are no NAs
    # --> expression levels of 0 obtain the lowest rank
    # turn into NAs to avoid detecting "lonely" interactions w/ only l or r
    if(!NA %in% rankLigs_curr | !NA %in% rankRecs_curr){
      print("Turnin into NAs")
      min_rank_lig <- min(rankLigs_curr)
      min_rank_rec <- min(rankRecs_curr)
      rankLigs_curr[rankLigs_curr == min_rank_lig] <- NA
      rankRecs_curr[rankRecs_curr == min_rank_rec] <- NA
    }
    
    # sums of NAs become NAs, as well -> "lonely" interactions become NAs
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


# wrapper for testing conditions during CCI calculation 


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## WRAPPER: CALCULATE SCORES
## wrapper function to calculate entire interactome in one step
## summarizes all steps from /04_cci_prep/02_cci_calculation
## all original functions are called from cci_functions_preparation.R 
## all alternative functions are defined above
## similar to calculate scores in cci_functions_perm

# returns CCI objects directly from a SCE object and lrdb

calculate_scores_wrap <- function(sce, lrdb, 
                                  cutoff = 0, 
                                  top_level, 
                                  assay_use,
                                  condition){
  
  # purpose: calculating final CCI object from pre-processed SCE object 
  
  sce <- sce              # SCE object, should be pre-processed
  lrdb <- lrdb            # ligand receptor database
  cutoff <- cutoff        # cutoff value between 0 and 1 (e.g. 0.2 = 20%)
  top_level <- top_level  # for leveling of ranks
  assay_use <- assay_use  # assay to use, "downsampled" by default
  condition <- condition  # type of pipeline (main, raw, no level, no cutoff)
  
  # extract species and age
  species <- sce$Species_ID[1]
  age <- sce$Age_ID[1]
  
  set.seed(37)
  
  #### STEP 1: GET PERCENTAGE OF EXPRESSING CELLS
  exp_gene_perc_df <- perc_expr_cells(sce, assay_use = assay_use)
  
  #### STEP 2: EXTRACT INTERACTION MATRIX
    # prepare extraction
  prep_list <- prepare_extraction(sce = sce, assay_use = assay_use)
  
  idents_list <- as.list(prep_list$idents)
  
    # in case of cutoff = 0, expr_genes_list will be all genes
  expr_genes_list <- lapply(X = idents_list, 
                            perc_df = exp_gene_perc_df, 
                            cutoff = cutoff, 
                            FUN = nexpr_gene_cutoff)
  
    # use prepared data to extract interaction matrix
    # this function was adjusted from Adrien's function interactionmatrix()
  
  names(expr_genes_list) <- prep_list$idents
  interaction_mat <- extract_matrix(counts = prep_list$countsmatrix,
                                    datasheet = prep_list$datasheet,
                                    expr_genes = expr_genes_list,
                                    lrdb = lrdb)
  
  #### STEP 3: INTERACTION SCORING
  datasheet <- prep_list$datasheet
  
  # this function was adjusted from Adrien's interactionranking() function
  interaction_score_list <- interaction_scoring_alt(
    int_df = interaction_mat, 
    datasheet = datasheet,
    top_level = top_level)

  #### STEP 4: MAKE INTERACTION LIST 
  interaction_list <- make_interaction_list(
    interaction_scores = interaction_score_list,
    lrdb = lrdb)
  
  interaction_list$Identities$species <- rep(species, 
                                             nrow(interaction_list$Identities))
  interaction_list$Identities$age <- rep(age, 
                                         nrow(interaction_list$Identities))
  interaction_list$Identities$condition <- rep(condition, 
                                               nrow(interaction_list$Identities))
  
  return(interaction_list)
}

