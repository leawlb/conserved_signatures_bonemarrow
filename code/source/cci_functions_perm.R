#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# functions for permutation testing as part of interactome analysis


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## PERM FUNCTION: STEP 0
## wrapper function to calculate entire interactome (cci object) in one step
## summarizes almost all steps from /04_cci_prep/02_cci_calculation except
## make_interaction_list() 
## for details please refer to /04_cci_prep/02_cci_calculation
## all functions are called from cci_functions_calculation.R
## required for permute_age()
## CALCULATE SCORES

# returns a data.table similar to CCI Score objects (data.table saves space)

calculate_scores <- function(
    sce, 
    min_perc, 
    lrdb,
    top_level){
  
  # calculate CCI-like object (interaction_list) from SCE object and LRDB
  
  sce <- sce              # preprocessed/downsampled SCE object
  min_perc <- min_perc    # minimum % of cells required to express a l/r gene 
  lrdb <- lrdb            # ligand receptor database
  top_level <- top_level  # top rank for better comparability between ranks
  
  species <- sce$Species_ID[1]
  age <- sce$Age_ID[1]
  
  print(assays(sce))
  
  set.seed(37)
  
  #### STEP 1: GET PERCENTAGE OF EXPRESSING CELLS
  exp_gene_perc_df <- perc_expr_cells(sce, assay_use = "downsampled")

  #### STEP 2: EXTRACT INTERACTION MATRIX
    # prepare extraction
  prep_list <- prepare_extraction(sce = sce, assay_use = "downsampled")
    # get identities (= cell types)
  idents_list <- as.list(prep_list$idents)

    # get list of genes expressed in at least min_perc % of cells per identity
  expr_genes_list <- lapply(X = idents_list, 
                            perc_df = exp_gene_perc_df, 
                            cutoff = min_perc, 
                            FUN = nexpr_gene_cutoff)
  
    # use prepared data to extract interaction matrix
    # this function was adjusted from Adrien's function interactionmatrix()
  names(expr_genes_list) <- prep_list$idents
  int_df <- extract_matrix(counts = prep_list$countsmatrix,
                           datasheet = prep_list$datasheet,
                           expr_genes = expr_genes_list,
                           lrdb = lrdb)
  
  #print(table(is.na(int_df)))
  
  #### STEP 3: INTERACTION SCORING
  datasheet <- prep_list$datasheet
  
    # this function was adjusted from Adrien's interactionranking() function
  interaction_score_list <- interaction_scoring(
    int_df = int_df, 
    datasheet = datasheet,
    top_level = top_level)
  
  # trying to save space
  remove(idents_list)
  remove(expr_genes_list)
  remove(datasheet)
  remove(int_df)
  remove(top_level)
  remove(prep_list)
  remove(min_perc)
  remove(sce)
  
  #print(table(is.na(interaction_score_list$Score)))
  
  #### STEP 4: MAKE INTERACTION LIST 
  #interaction_list <- make_interaction_list(
  #  interaction_scores = interaction_score_list,
  #  lrdb = lrdb)
    # not required, omitted to save space
  
  score_dt <- data.table(interaction_score_list$Score)
  remove(interaction_score_list)
  cci_like <- list("Score" = score_dt, "Identities" = c("species" = species,
                                                        "age" = age)) 

  return(cci_like)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


## AGE PERMUTATION 


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## PERM FUNCTION: STEP 1
## permute age labels and calculate scores for "young" and "old" SCE objects
## PERMUTE LABELS: AGE

# returns two permuted CCI objects, one for each age
# use with mclapply on list of iterations --> list of CCIs
# requires calculate_scores()

permute_age <- function(
    iterations_vec,
    sce_old, 
    sce_yng,
    lrdb,
    min_perc,
    top_level){
  
  # purpose: permute age labels and calculate CCI for perm yng/old SCE objects
  
  i <- iterations_vec      # vector from 1 to nr of iterations
  sce_old <- sce_old       # old SCE object for perm testing
  sce_yng <- sce_yng       # yng SCE object for perm testing
  lrdb <- lrdb             # ligand receptor database
  min_perc <- min_perc     # minimum % of cells required to express a l/r gene 
  top_level <- top_level   # top rank for better comparability between ranks
  
  species <- sce_old$Species_ID[1]
  stopifnot(sce_old$Species_ID[1] == sce_yng$Species_ID[1])
  
  #### PERMUTE LABELS 
  sce_perm <- cbind(sce_old, sce_yng)
  
  sce_perm$age_permuted <- vector(length = ncol(sce_perm))
  
    # permute age labels within cell types (identities), otherwise ct labels 
    # will also be permuted
  cts <- levels(sce_perm$Identity)
  for(c in 1:length(cts)){
    sce_ct <- sce_perm[,sce_perm$Identity == cts[c]]
    
    # choose seed depending on iteration and ct
    set.seed(i*7144+c+11)
    perm_ages <- sample(c(1:ncol(sce_ct)), size = ncol(sce_ct))
    
    set.seed(37) # reset seed to basic value
    age_permuted <- sce_ct$Age_ID[perm_ages]
    sce_perm[,sce_perm$Identity == cts[c]]$age_permuted <- age_permuted
  
  }
  print(table(sce_perm$Age_ID, sce_perm$age_permuted))
  
  sce_perm_yng <- sce_perm[,sce_perm$age_permuted == "yng"]
  sce_perm_old <- sce_perm[,sce_perm$age_permuted == "old"]

  #### CALCULATE CCI OBJECTS OF PERMUTED SCE OBJECTS WITH calculate_scores()
    # separately for "yng" and "old" perm SCEs
  cci_like_yng <- calculate_scores(
    sce = sce_perm_yng, 
    lrdb = lrdb,
    min_perc = min_perc,
    top_level = top_level)

  cci_like_old <- calculate_scores(
    sce = sce_perm_old, 
    lrdb = lrdb,
    min_perc = min_perc,
    top_level = top_level)    

  return_list <- list("yng" = cci_like_yng, "old" = cci_like_old)
  return(return_list)
  
  # for simplicity, will still be referred to as "cci"
  
}

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## PERM FUNCTION: STEP 2
## Convert lists of permuted CCIs per age into list of permuted score df:
## CONVERT_PERM_DF()

# returns a list of dataframes with rows = interactions and cols = iterations
# containing the permuted scores
# one list item per ctp
# use with apply for lists from different ages 

convert_perm_df <- function(cci_list_age, lrdb){
  
  # purpose: convert lists of permuted CCIs into dfs with permuted scores
  
  cci_list_age <- cci_list_age  # list permuted cci object per age
  lrdb <- lrdb                  # ligand receptor database
  
  # obtain interactions, always same sequence as lrdb$lr_pair
  interactions <- lrdb$lr_pair
  # obtain cell type pairs, always the same ctps for each permuted CCI pair
  ctps <- colnames(cci_list_age[[1]]$Score)
  
  # for each cell type pair ctp
  score_df_list <- lapply(as.list(ctps), 
                          function(ctp = as.list(ctps),
                                   interactions = interactions, 
                                   cci_list_age = cci_list_age){
    
    ctp <- ctp
    interactions <- interactions # vector of interactions as rownames for df
    cci_list_age <- cci_list_age # list of ccis

    # for each permuted cci object in cci_list_age, cci
    temp_df_list <- lapply(cci_list_age, 
                           function(cci = cci_list_age, 
                                    interactions = interactions, 
                                    ctp = ctp){
      
      cci <- cci                     # current permuted CCI object
      interactions <- interactions   # vector of interactions as rownames for df
      ctp <- ctp                     # current cell type pair

      # convert to data.frame for better handling
      score <- as.data.frame(cci$Score)
      colnames(score) <- colnames(cci$Score)
      
      # df with rows = interactions, fill with perm scores of current ctp 
      temp_df <- data.frame(row.names = interactions)
      temp_df[,1] <- score[,colnames(score) == ctp]
  
      return(temp_df)
      
    }, interactions, ctp)

    names(temp_df_list) <- names(cci_list_age)
    
    suppressMessages(score_df <- bind_cols(temp_df_list))
    colnames(score_df) <- names(cci_list_age)
    # this is now one score df for one ctp
    return(score_df)
    
  }, interactions, cci_list_age)
  
  # a list of score dfs of length = n ctps
  names(score_df_list) <- ctps
  
  # return for current age
  return(score_df_list)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## PERM FUNCTION: STEP 2
## this function calculates important statistics after permutation of age labels
## it obtains score deltas (score_old - score_yng), and permuted score deltas. 
## It estimates p values, adjusts them, and evaluates normal distribution of 
## permuted score deltas obtained from permute_age()
## STATS_PERM_AGE_PREP

# required for stats_perm_age()
# returns per cell type pair:
# a) one row of a df containing various stats with rows = interactions, cols = stats
# b) the permutated score deltas (sds) for visualisation

# use lapply on list of interactions within one ctp (cell type pair)
# 1. -> list of one-row dfs, one df per interaction
# bind_rows, then use lapply on list of ctps 
# 2. -> list of dfs/sds per ctp

stats_perm_age_prep <- function(
    int, 
    ctp, 
    iterations, 
    cci_yng,
    cci_old,
    score_dft_yng, 
    score_dft_old,
    min_sp_age_perc){
  
  # purpose: calculate stats from score deltas obtained from permutation per ctp
  
  int <- int                          # ligand receptor interaction pair
  ctp <- ctp                          # current cell type pair
  iterations <- iterations            # nr of iterations used in previous perm
  cci_yng <- cci_yng                  # young CCI object (non-permuted)
  cci_old <- cci_old                  # old CCI object (non-permuted)
  score_dft_yng <- score_dft_yng      # perm sds, from "yng" perm objects
  score_dft_old <- score_dft_old      # perm sds, from "old" perm objects
  min_sp_age_perc <- min_sp_age_perc  # min % of successful perms (sp)
  # most CCI scores are NAs (= not detected/not "successful")
  # for some ints, a low number of "successfully" calculated scores (<100) 
  # is possible, but this is not enough to calculate p vals
  
  all_ints <- colnames(score_dft_yng)
  
  #### PREPARATION
    # empty res_df to fill
    # only one row is used in this function, corresponding to current int
    # a list of one-row dfs will later be combined by bind_rows()
  res_df_empty <- data.frame(
    row.names = all_ints,
    "score_yng" = vector(length = length(all_ints)),
    "score_old" = vector(length = length(all_ints)),
    "shared" = vector(length = length(all_ints)),
    "score_delta" = vector(length = length(all_ints)),
    
    # pval 
    "p_val" = vector(length = length(all_ints)),
    "p_adj" = vector(length = length(all_ints)),
    
    # perm score delta stats
    "perm_scores_delta_mean" = vector(length = length(all_ints)),
    "perm_scores_delta_median" = vector(length = length(all_ints)),
    "perm_scores_delta_length" = vector(length = length(all_ints)),
    "perm_scores_delta_sd" = vector(length = length(all_ints)),
    "perm_scores_delta_iqr" = vector(length = length(all_ints)),
    "perm_scores_delta_min" = vector(length = length(all_ints)),
    "perm_scores_delta_max" = vector(length = length(all_ints)),
    "perm_scores_delta_nr_obs" = vector(length = length(all_ints)),
    "perm_scores_delta_min_centered" = vector(length = length(all_ints)),
    "perm_scores_delta_max_centered" = vector(length = length(all_ints)),
    "score_delta_centered" = vector(length = length(all_ints)),
    
    # norm distribution test 
    "p_val_shapiro" = vector(length = length(all_ints)),
    "p_adj_shapiro" = vector(length = length(all_ints)))
  
  i <- which(all_ints == int)
  temp_df <- res_df_empty[i,]
  
    # get the vector of permuted scores for current interaction nr i
  perm_scores_yng <- score_dft_yng[,i]
  perm_scores_old <- score_dft_old[,i]

    # get the scores from the cci object
  pos_ctp_yng <- which(colnames(cci_yng$Score) == ctp)
  score_yng <- cci_yng$Score[which(rownames(cci_yng$Score) == int), pos_ctp_yng]
  temp_df$score_yng[1] <- score_yng
  
  pos_ctp_old <- grep(ctp, colnames(cci_old$Score))
  score_old <- cci_old$Score[which(rownames(cci_old$Score) == int), pos_ctp_old]
  temp_df$score_old[1] <- score_old
  
  #### GET SCORE DELTAS AND PERMUTED SCORE DELTAS
    # only for permutations where both young and old scores are not NAs
    # and fill info into df 
  if(!is.na(score_old) & !is.na(score_yng)){

    temp_df$shared[1] <- TRUE
    
    score_delta <- score_old - score_yng
    temp_df$score_delta[1] <- score_delta
    
    perm_score_deltas <- perm_scores_old - perm_scores_yng
    
    temp_df$perm_scores_delta_length[1] <- length(
      which(!is.na(perm_score_deltas)))
    
      # add stats only if at least a certain nr of permutations was successful
      # this prevents calculating means etc. for NAs or very low nr of sds only
    min_sp_age <- iterations * min_sp_age_perc
    if(length(which(!is.na(perm_score_deltas))) >= min_sp_age){
      
      temp_df$perm_scores_delta_mean[1] <- mean(
        perm_score_deltas[!is.na(perm_score_deltas)])
      temp_df$perm_scores_delta_median[1] <- median(
        perm_score_deltas[!is.na(perm_score_deltas)])
      temp_df$perm_scores_delta_sd[1] <- sd(
        perm_score_deltas[!is.na(perm_score_deltas)])
      temp_df$perm_scores_delta_iqr[1] <- IQR(
        perm_score_deltas[!is.na(perm_score_deltas)])
      temp_df$perm_scores_delta_min[1] <- min(
        perm_score_deltas[!is.na(perm_score_deltas)])
      temp_df$perm_scores_delta_max[1] <- max(
        perm_score_deltas[!is.na(perm_score_deltas)])
    
      # get pval, which is the number of permuted deltas that are as extreme or
      # more extreme than the test statistic = score_delta
      # test using difference from median, because median might not be exactly 0

      # remove NAs from calculation
      perm_score_deltas_nna <- perm_score_deltas[!is.na(perm_score_deltas)]
      
      ## TODO: check if this is correct!
      # center perm distribution and test statistic around 0 
      perm_sd_centered <- perm_score_deltas_nna - median(perm_score_deltas_nna)
      score_delta_centered <- score_delta - median(perm_score_deltas_nna)
      
      # obtain nr of observed permuted values that is at least as or more  
      # extreme as test statistic score delta
      nr_obs <- length(
        which(abs(perm_sd_centered) >= abs(score_delta_centered)))
      
      temp_df$perm_scores_delta_nr_obs[1] <- nr_obs
      temp_df$perm_scores_delta_min_centered[1] <- min(perm_sd_centered)
      temp_df$perm_scores_delta_max_centered[1] <- max(perm_sd_centered)
      temp_df$score_delta_centered[1] <- score_delta_centered
      
      pval <- nr_obs/iterations
      temp_df$p_val[1] <- pval
    
      if(pval == 0){
        warning(paste0("Pval is exactly 0 for ", int))
      }
    
    }else{
  
      # set to NA if nr of successful perm is (too) low
      temp_df$perm_scores_delta_mean[1] <- NA
      temp_df$perm_scores_delta_median[1] <- NA
      temp_df$perm_scores_delta_sd[1] <- NA
      temp_df$perm_scores_delta_iqr[1] <- NA
      temp_df$perm_scores_delta_min[1] <- NA
      temp_df$perm_scores_delta_max[1] <- NA
      temp_df$p_val[1] <- NA 
      
    }
    
    # test normal distribution of permuted score deltas
    if(length(unique(perm_score_deltas)) > 1 & length(which(!is.na(perm_score_deltas))) > 3){
      temp_df$p_val_shapiro[1] <- shapiro.test(perm_score_deltas)[2]$p.value
    }else{
      temp_df$p_val_shapiro[1] <- NA
    }
  }else{
    # if one or both of the actual scores is NA, no delta can be calculated
    # fill df with NAs
    temp_df$shared[1] <- FALSE
    temp_df$score_delta[1] <- NA
    temp_df$p_val[1] <- NA
    temp_df$p_adj[1] <- NA
    temp_df$perm_scores_delta_mean[1] <- NA
    temp_df$perm_scores_delta_median[1] <- NA
    temp_df$perm_scores_delta_length[1] <- NA
    temp_df$perm_scores_delta_sd[1]  <- NA
    temp_df$perm_scores_delta_iqr[1]  <- NA
    temp_df$perm_scores_delta_min[1]  <- NA
    temp_df$perm_scores_delta_max[1]  <- NA
    
    temp_df$perm_scores_delta_nr_obs[1] <- NA
    temp_df$perm_scores_delta_min_centered[1] <- NA
    temp_df$perm_scores_delta_max_centered[1] <- NA
    temp_df$score_delta_centered[1] <- NA
    
    temp_df$p_val_shapiro[1] <- NA
    temp_df$p_adj_shapiro[1] <- NA
    
    perm_score_deltas <- rep(NA, iterations)
  }
  
  return_list <- list(temp_df, perm_score_deltas)
  names(return_list) <- c("df", "perm_score_deltas")
  return(return_list)
  
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## PERM FUNCTION: STEP 2
## extract the stats from age-permuted CCI objects per cell type pair
## extract age-permuted score deltas for visualisation
## STATS_PERM_AGE

# this function returns a list of dfs and permuted score delta vectors per
# cell type pair (ctp)
# requires stats_perm_age_prep()

stats_perm_age <- function(
    ctp,
    perm_score_list_yng, 
    perm_score_list_old,
    cci_yng,
    cci_old,
    iterations,
    min_sp_age_perc){
  
  ctp <- ctp                                  # current cell type pair
  perm_score_list_yng <- perm_score_list_yng  # young list from perm_age()
  perm_score_list_old <- perm_score_list_old  # old list from perm_age()
  cci_yng <- cci_yng                          # young CCI object (non-perm)
  cci_old <- cci_old                          # old CCI object (non-perm)
  iterations <- iterations                    # nr of iterations previous perm
  min_sp_age_perc <- min_sp_age_perc          # min % of successful iterations 
  # most CCI scores are NAs (= not detected/not "successful")
  # for some ints, a low number of "successfully" calculated scores (<100) 
  # is possible, but this is not enough to calculate p vals
  
  species <- cci_yng$Species[1]
  
  # transpose score_dfs for stats_perm_age_prep()
  score_dft_yng <- as.data.frame(t(perm_score_list_yng[[ctp]]$df))
  score_dft_old <- as.data.frame(t(perm_score_list_old[[ctp]]$df))
  
  int_list <- as.list(colnames(score_dft_yng))
  
  print("starting pval estimation")
  
  res_list <- lapply(
    X = int_list,
    ctp = ctp,
    iterations = iterations, 
    min_sp_age_perc = min_sp_age_perc,
    cci_yng = cci_yng,
    cci_old = cci_old,
    score_dft_yng = score_dft_yng, 
    score_dft_old = score_dft_old,
    FUN = stats_perm_age_prep
  )
  
  names(res_list) <- unlist(int_list)
  
  # separate res_list into stat_df slices (one-row dfs) and perm sd slices
  res_df_slices <- list()
  perm_score_delta_slices <- list()
  for(i in unlist(int_list)){
    res_df_slices[[i]] <- res_list[[i]]$df
    perm_score_delta_slices[[i]] <- res_list[[i]]$perm_score_deltas
  }

  print("finished pval estimation")
  
  # combine res_df slices into one stats_df per current ctp
  stats_df <- bind_rows(res_df_slices)
  
  # adjust p value within df
  stats_df$p_adj <- p.adjust(stats_df$p_val, method = "fdr")
  stats_df$p_adj_shapiro <- p.adjust(stats_df$p_val_shapiro, method = "fdr")
  
  # add metadata
  stats_df$species <- rep(species, nrow(stats_df))
  stats_df$ctp <- rep(ctp, nrow(stats_df))
  
  # finish perm_score_delta df for visualisation
  vis_df <- data.frame(iterations =  paste0("iteration_", c(1:iterations)))
  for(i in 1:length(perm_score_delta_slices)){
    vis_df[,1+i] <- perm_score_delta_slices[[i]]
    colnames(vis_df)[1+i] <- names(perm_score_delta_slices)[i]
  }
  vis_df <- column_to_rownames(vis_df, var = "iterations")
  
  return_list <- list(stats_df, vis_df)
  names(return_list) <- c("res_df", "vis_df")
  return(return_list)
}


#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


## CELL TYPE PERMUTATION


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## PERMUTATION OF CELL TYPE LABELS

permute_cts <- function(iterations_vec, sce_cond,
                        lrdb, min_perc, top_level){
  
  # sce to permute
  sce_perm <- sce_cond
  i <- iterations_vec
  #print(i)
  
  # permute labels using a seed dependent on nr of iteration
  # for emitters and receivers separately
  # the number of cells per cell type will always be the same 
  
  set.seed((i*317))
  # get permutation vectors for each Assignment separately
  perm_emis <- sample(c(1:length(which(sce_perm$Assignment == "emitter"))), 
                      size = length(which(sce_perm$Assignment == "emitter")))
  perm_recs <- sample(c(1:length(which(sce_perm$Assignment == "receiver"))),
                      size = length(which(sce_perm$Assignment == "receiver")))
  
  # record status to compare
  sce_perm$assignment_before <- sce_perm$Assignment
  sce_perm$identity_before <- sce_perm$Identity
  
  # permute
  sce_perm$Identity[
    which(sce_perm$Assignment == "emitter")] <- sce_perm$Identity[
      which(sce_perm$Assignment == "emitter")][perm_emis]
  sce_perm$Identity[
    which(sce_perm$Assignment == "receiver")] <- sce_perm$Identity[
      which(sce_perm$Assignment == "receiver")][perm_recs]
  
  #print(table(sce_perm$Assignment, sce_perm$assignment_before))
  #print(table(sce_perm$Identity, sce_perm$identity_before))
  #print(table(sce_perm$Identity))
  
  # quickly calculate CCI in one step 
  int_list <- calculate_scores(sce_perm, 
                               lrdb = lrdb,
                               min_perc = min_perc,
                               top_level = top_level, 
                               age = sce_perm$Age_ID[1])
  
  return(int_list)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## get important/interesting stats for ct permutations

stats_ints_cts <- function(int, ctp, iterations, min_sp_cts,
                           cci_cond,
                           score_dft_cond){
  
  # to iterate through interactions
  all_ints <- colnames(score_dft_cond)
  
  # empty res_df to fill, but only one row will be used in this function
  res_df_empty <- data.frame(
    row.names = all_ints,
    "score" = vector(length = length(all_ints)),

    # pval 
    "p_val" = vector(length = length(all_ints)),
    "p_adj" = vector(length = length(all_ints)),
    
    # perm score stats
    "perm_scores_mean" = vector(length = length(all_ints)),
    "perm_scores_median" = vector(length = length(all_ints)),
    "perm_scores_length" = vector(length = length(all_ints)),
    "perm_scores_sd" = vector(length = length(all_ints)),
    "perm_scores_iqr" = vector(length = length(all_ints)),
    "perm_scores_min" = vector(length = length(all_ints)),
    "perm_scores_max" = vector(length = length(all_ints)),
    "perm_scores_nr_obs" = vector(length = length(all_ints)),
    "perm_scores_min_centered" = vector(length = length(all_ints)),
    "perm_scores_max_centered" = vector(length = length(all_ints)),
    "score_centered" = vector(length = length(all_ints)),
    
    # norm distribution test
    "p_val_shapiro" = vector(length = length(all_ints)),
    "p_adj_shapiro" = vector(length = length(all_ints)))
  
  i <- which(all_ints == int)
  temp_df <- res_df_empty[i,]
  
  # get the vector of permuted scores of interaction i
  perm_scores <- score_dft_cond[,i]

  # get the score from the cci object at position i, ctp
  pos_ctp <- which(colnames(cci_cond$Score) == ctp)
  score <- cci_cond$Score[which(rownames(cci_cond$Score) == int), pos_ctp]
  temp_df$score[1] <- score
  
  # calculate pval only if score is not NA
  if(!is.na(score)){
    
    # only if certain percentage of permutations was successful
    # also prevents calculating means etc. for rare all NAs
    # min_sp_cts can be a little smaller than  as cell types are 
    # more different from each other than same cell types at different ages
    min_sp_cts <- iterations * min_sp_ct_perc
    if(length(which(!is.na(perm_scores))) >= min_sp_cts){
      
      temp_df$perm_scores_mean[1] <- mean(perm_scores[!is.na(perm_scores)])
      temp_df$perm_scores_median[1]  <- median(perm_scores[!is.na(perm_scores)])
      temp_df$perm_scores_length[1]  <- length(which(!is.na(perm_scores)))
      temp_df$perm_scores_sd[1]  <- sd(perm_scores[!is.na(perm_scores)])
      temp_df$perm_scores_iqr[1]  <- IQR(perm_scores[!is.na(perm_scores)])
      temp_df$perm_scores_min[1]  <- min(perm_scores[!is.na(perm_scores)])
      temp_df$perm_scores_max[1]  <- max(perm_scores[!is.na(perm_scores)])
    
      # get pval, which is the number of permuted deltas that are as extreme or
      # more extreme than the test statistic = score_delta

      # remove NAs from calculation
      perm_scores_tc <- perm_scores[!is.na(perm_scores)]
      
      # center distribution and values around 0 
      perm_scores_centered <- perm_scores_tc - median(perm_scores_tc)
      score_centered <- score - median(perm_scores_tc)
      
      # obtain nr of values that is at least as or more extreme than test statistic
      nr_obs <- length(which(abs(perm_scores_centered) >= 
                               abs(score_centered)))
      
      # record
      temp_df$perm_scores_nr_obs[1] <- nr_obs
      temp_df$perm_scores_min_centered[1] <- min(perm_scores_centered)
      temp_df$perm_scores_max_centered[1] <- max(perm_scores_centered)
      temp_df$score_centered[1] <- score_centered
      
      pval <- nr_obs/iterations
      temp_df$p_val[1] <- pval
      
      if(pval == 0){
        warning(paste0("Pval is exactly 0 for ", int))
      }
    }else{
      temp_df$p_val[1] <- NA
    }
    
    # test for normal distribution of permuted score deltas
    # only if score is not NA
    if(length(unique(perm_scores)) > 3){
      temp_df$p_val_shapiro[1] <- shapiro.test(perm_scores)[2]$p.value
    }else{
      temp_df$p_val_shapiro[1] <- NA
    }
  }else{
    temp_df$score[1] <- NA
    temp_df$p_val[1] <- NA
    temp_df$p_adj[1] <- NA
    temp_df$perm_scores_mean[1] <- NA
    temp_df$perm_scores_median[1] <- NA
    temp_df$perm_scores_length[1] <- NA
    
    temp_df$perm_scores_sd[1]  <- NA
    temp_df$perm_scores_iqr[1]  <- NA
    temp_df$perm_scores_min[1]  <- NA
    temp_df$perm_scores_max[1]  <- NA
    
    temp_df$perm_scores_nr_obs[1]  <- NA
    temp_df$perm_scores_min_centered[1]  <- NA
    temp_df$perm_scores_max_centered[1]  <- NA
    temp_df$score_centered[1]  <- NA
    
    temp_df$p_val_shapiro[1] <- NA
    temp_df$p_adj_shapiro[1] <- NA
    
    perm_scores <- rep(NA, iterations)
  }
  
  #print(table(is.na(perm_scores)))
  return_list <- list(temp_df, perm_scores)
  names(return_list) <- c("df", "perm_scores")
  return(return_list)
  
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

# obtain a list of interesting/important stats for ct permutation
# requires stats_ints_cts 

get_interaction_stats_cts <- function(ctp, species, age, min_sp_cts,
                                      perm_score_list_cond, 
                                      cci_cond,
                                      iterations){
  
  # transpose score_dfs
  score_dft_cond <- as.data.frame(t(perm_score_list_cond[[ctp]]$df))

  int_list <- as.list(colnames(score_dft_cond))
  
  print("starting pval estimation")
  
  print(ctp)
  
  res_list <- lapply(X = int_list, ctp = ctp, iterations = iterations, 
                     min_sp_cts = min_sp_cts,
                     cci_cond = cci_cond, 
                     score_dft_cond = score_dft_cond, 
                     FUN = stats_ints_cts)
  names(res_list) <- unlist(int_list)
  
  res_df_slices <- list()
  perm_score_slices <- list()
  for(i in unlist(int_list)){
    res_df_slices[[i]] <- res_list[[i]]$df
    perm_score_slices[[i]] <- res_list[[i]]$perm_scores
  }
  #print(names(res_list[[1]]))
  
  print("finished pval estimation")
  
  # finish res_df containing pvals
  res_df <- bind_rows(res_df_slices)
  
  res_df$p_adj <- p.adjust(res_df$p_val, method = "fdr")
  res_df$p_adj_shapiro <- p.adjust(res_df$p_val_shapiro, method = "fdr")
  
  res_df$species <- rep(species, nrow(res_df))
  res_df$age <- rep(age, nrow(res_df))
  res_df$ctp <- rep(ctp, nrow(res_df))
  
  # finish perm_scores df for visualisation
  vis_df <- data.frame(iterations =  paste0("iteration_", c(1:iterations)))
  for(i in 1:length(perm_score_slices)){
    vis_df[,1+i] <- perm_score_slices[[i]]
    colnames(vis_df)[1+i] <- names(perm_score_slices)[i]
  }
  vis_df <- column_to_rownames(vis_df, var = "iterations")
  
  return_list <- list(res_df, vis_df)
  names(return_list) <- c("res_df", "vis_df")
  return(return_list)
}
