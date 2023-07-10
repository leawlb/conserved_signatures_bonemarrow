# functions for analysis of cell type interactomes
# for detailed explanations see .Rmd files for cpi analysis.

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# METRICS FUNCTIONS 
# required to calculate important metrics of cpi objects
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

extract_ctp_info <- function(cpi){
  # cell type pair interactome object
  
  # empty objects
  return_list <- list()
  return_df <- data.frame(
    "ct_pairs" = colnames(cpi$Score)
  )
  return_df$nr_ints <- vector(length = nrow(return_df))
  return_df$emitter <- vector(length = nrow(return_df))
  return_df$receiver <- vector(length = nrow(return_df))
  
  # for each ctp in the cpi object, perform separately
  # this will yield ncol dataframes per cpi object in a large return_list
  for(i in 1:ncol(cpi$Score)){
    temp_df <- data.frame(
      "ct_pair" = rep(colnames(cpi$Score)[i], nrow(cpi$Score)),
      "interaction" = rownames(cpi$Score),
      "score" = cpi$Score[,i],
      "interaction_type" = cpi$Interactions$interaction_type,
      "emitter" = rep(cpi$Celltypes$emitter[i], nrow(cpi$Score)),
      "ligandrank" = cpi$Ligandrank[,i],
      "receiver" = rep(cpi$Celltypes$receiver[i], nrow(cpi$Score)),
      "receptorrank" = cpi$Receptorrank[,i]
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
  names(return_list) <- colnames(cpi$Score)
  return_df <- return_df[order(return_df$nr_ints, decreasing = TRUE),]
  
  # the resulting object contains an overview df of all ctp interactions, 
  # then separate dfs for each separate ctp
  return(list("overview" = return_df, 
              "separate" = return_list))
}

# this function extracts important info on interactions at CELL TYPE PAIR level
# can only be used if cutoff = TRUE and level = TRUE and no subsetting required
# otherwise, use extract_ctp_info_raw()

#-------------------------------------------------------------------------------

extract_ctp_info_raw <- function(cpi, 
                                 ct_pairs = NULL, 
                                 itypes = NULL, 
                                 cutoff = FALSE, 
                                 level = FALSE){    
  # cpi = cell type pair interactome object
  # ct_pairs = vector of ctp names to keep (subset)
  # itypes = vector of int types to keep (subset)
  # cutoff = indicator whether gene cut off
  # level = indicator whether leveling of ranks was implemented
  
  #---------subset the cpi according to ct_pairs and interaction types-----------                                 
  if(is.null(ct_pairs) == TRUE){  
    ct_pairs <- colnames(cpi$Score)[1:ncol(cpi$Score)]
  }else if(is.null(ct_pairs) == FALSE){
    cpi$Score <- cpi$Score[,!is.na(match(colnames(cpi$Score), ct_pairs))]
    cpi$Ligandrank <- cpi$Ligandrank[
      ,!is.na(match(colnames(cpi$Ligandrank), ct_pairs))]
    cpi$Receptorrank <- cpi$Receptorrank[
      ,!is.na(match(colnames(cpi$Receptorrank), ct_pairs))]
  }
  if(is.null(itypes) == FALSE){
    keep_rows <- vector(mode = "numeric", length = 0)
    for(i in length(itypes)){
      keep_rows <- c(which(cpi$Interactions$interaction_type == itypes[i]))
    }
    cpi$Score <- cpi$Score[keep_rows,]
  }
  #----------------------------------------------------------------------------- 
  
  # keep lri names when they can be counted as "detected"
  # "detected" is defined differently for each object depending on conditions 
  
  if(cutoff == FALSE && level == FALSE){  # raw 
    # lris are counted as detected if both l and r ranks are larger than 1
    
    temp_df <- cpi$Score # to borrow nrow and ncol
    temp_list_ints <- list()
    # iterate through cols AND rows to make sure both l and r are expressed
    for(j in 1:ncol(cpi$Score)){
      for(i in 1:nrow(cpi$Score)){
        if(cpi$Receptorrank[as.numeric(i), as.numeric(j)] > 1  &&
           cpi$Ligandrank[as.numeric(i), as.numeric(j)] > 1){
          temp_df[i,j] <- rownames(cpi$Receptorrank)[i]
        }else{
          temp_df[i,j] <- NA
        }
      }
      # now for each ctp, put only detected interactions into list
      temp_list_ints[[j]] <- temp_df[!is.na(temp_df[,j]), j]
    }
    names(temp_list_ints) <- colnames(cpi$Score)
    
  }else if(cutoff == FALSE){  # nct
    # lris are counted as detected if score is larger than 0
    
    temp_list_ints <- list()
    for(j in 1:ncol(cpi$Score)){
      temp_list_ints[[j]] <- rownames(cpi$Score)[which(cpi$Score[,j] > 0)]
    }
    names(temp_list_ints) <- colnames(cpi$Score)
    
  }else if(cutoff == TRUE){ # FULL, nlv, nrm
    # lris are counted as detected if score is not NA
    
    temp_list_ints <- list()
    for(j in 1:ncol(cpi$Score)){
      temp_list_ints[[j]] <- rownames(cpi$Score)[which(!is.na(cpi$Score[,j]))]
    }
  }
  
  temp_df1_list <- list()
  for(i in 1:length(temp_list_ints)){
    # make a df to return for each CTP
    temp_df1 <- data.frame(
      ct_pair = rep(colnames(cpi$Score)[i], length(temp_list_ints[[i]])),
      interaction = temp_list_ints[[i]],
      score = vector(length = length(temp_list_ints[[i]])),
      interaction_type = vector(length = length(temp_list_ints[[i]])),
      emitter =  vector(length = length(temp_list_ints[[i]])),
      ligandrank =  vector(length = length(temp_list_ints[[i]])),
      receiver =  vector(length = length(temp_list_ints[[i]])),
      receptorrank =  vector(length = length(temp_list_ints[[i]]))
    )
    
    for(j in 1:nrow(temp_df1)){
      temp_df1$score[j] <- cpi$Score[
        which(rownames(cpi$Score) %in% temp_df1$interaction[j]),i]
      temp_df1$interaction_type[j] <- cpi$Interaction$interaction_type[
        which(cpi$Interaction$interaction_pair %in% temp_df1$interaction[j])]
      temp_df1$ligandrank[j] <- cpi$Ligandrank[
        which(rownames(cpi$Ligandrank) %in% temp_df1$interaction[j]),i]
      temp_df1$receptorrank[j] <- cpi$Receptorrank[
        which(rownames(cpi$Receptorrank) %in% temp_df1$interaction[j]),i]
    }
    
    temp_df1$emitter <- gsub("[&][[:print:]]+", "", temp_df1$ct_pair)
    temp_df1$receiver <- gsub("[[:print:]]+[&]", "", temp_df1$ct_pair)
    
    temp_df1_list[[i]] <- temp_df1
  }
  names(temp_df1_list) <- colnames(cpi$Score)
  
  # add only number of interactions to a separate overview df
  temp_df2 <- data.frame(row.names = colnames(cpi$Score))
  for(i in 1:nrow(temp_df2)){
    temp_df2[i,1] <- length(temp_list_ints[[i]])
  }  
  colnames(temp_df2)[1] <- "nr_ints"
  temp_df2 <- rownames_to_column(temp_df2, var = "ct_pairs")
  temp_df2 <- temp_df2[order(temp_df2$nr_ints, decreasing = TRUE),]
  
  temp_df2$emitter <- gsub("&.*", "", temp_df2$ct_pairs)
  temp_df2$receiver <- gsub(".*&", "", temp_df2$ct_pairs)
  
  # add all lri to list
  return_list <- list(
    overview = temp_df2,
    separate = temp_df1_list
  )
  names(return_list$separate) <- ct_pairs
  
  return(return_list)

}  

# this function extracts important info on interactions at CELL TYPE PAIR level
# it can still be used for cpi objects calculated under different conditions
# and to subset the cpi objects by interaction types or cell type 
# this is still functional, but depreciated for cutoff = TRUE, level = TRUE cpi


#-------------------------------------------------------------------------------

extract_ct_info <- function(ctpints){
  
  # put all info from ctpints lists into one big dataframe
  ctpints_total <- ctpints[[2]][[1]]
  for(i in 2:length(ctpints[[2]])){

    ctpints_total <- rbind(ctpints_total, ctpints[[2]][[i]])
  }
  
  # then subset by cell types, first emitters
  nr_emitters <- length(unique(ctpints_total$emitter))
  nr_receivers <- length(unique(ctpints_total$receiver))
  
  ctpints_total_list<- list()
  for(i in 1:nr_emitters){
    ctpints_total_list[[i]] <- ctpints_total[
      ctpints_total$emitter == unique(ctpints_total$emitter)[i],]
    ctpints_total_list[[i]]$celltype <- unique(ctpints_total$emitter)[i]
    ctpints_total_list[[i]]$annotation <- "emitter"
  }
  
  # then receivers
  for(i in 1:nr_receivers){
    
    j <- nr_emitters + i
    ctpints_total_list[[j]] <- ctpints_total[
      ctpints_total$receiver == unique(ctpints_total$receiver)[i],]
    ctpints_total_list[[j]]$celltype <- unique(ctpints_total$receiver)[i]
    ctpints_total_list[[j]]$annotation <- "receiver"
  }
  
  # now use this data to make a ctints_list with all relevant information
  ctints_list <- lapply(ctpints_total_list, function(ctpints_total){
    
    ctints <- data.frame(
      interactions = unique(ctpints_total$interaction),
      celltype = rep(ctpints_total$celltype, 
                     length = length(unique(ctpints_total$interaction))),
      annotation = rep(ctpints_total$annotation,
                       length = length(unique(ctpints_total$interaction)))
    )
    
    # add the relevant rank
    if(ctints$annotation[1] == "emitter"){
      ctints$rank <-  ctpints_total$ligandrank[
        match(ctints$interactions, ctpints_total$interaction)]
    }else if(ctints$annotation[1] == "receiver"){
      ctints$rank <-  ctpints_total$receptorrank[
        match(ctints$interactions, ctpints_total$interaction)]
    }
    
    ctints <- ctints[order(ctints$rank, decreasing = TRUE),]
    return(ctints)
  })
  
  names(ctints_list) <- c(unique(ctpints_total$emitter), 
                          unique(ctpints_total$receiver))
  
  # make an overview df
  overview_df <- data.frame(
    celltypes =  names(ctints_list)
  )
  
  for(i in 1:length(ctints_list)){
    overview_df$nr_ints[i] <- nrow(ctints_list[[i]])
    overview_df$annotation[i] <- ctints_list[[i]]$annotation[1]
  }
  overview_df <- overview_df[order(overview_df$nr_ints, decreasing = TRUE),]
  
  return(list(overview = overview_df, 
              separate = ctints_list))
  
}


# this function uses the ctpints list to extract info on a cell type level

#-------------------------------------------------------------------------------

extract_lrs_info <- function(ctints){
  # object with info on interactions on CELL TYPE level from extract_ct_info()
  
  if("emitter" %in% ctints$annotation){
    # remove the receiver from the lri name to obtain only the ligand
    ctints$interactions <- gsub("[&][[:alnum:]]+", "", ctints$interactions)
    ligands <- unique(ctints$interactions)
    
    return_df <- data.frame(
      celltypes = rep(ctints$celltype[1], length(ligands)),
      annotations = rep("emitter", length(ligands)),
      sigmols = ligands
    )
  }
  if("receiver" %in% ctints$annotation){
    # remove the emitter from the lri name to obtain only the receptor
    ctints$interactions <- gsub("[[:alnum:]]+[&]", "", ctints$interactions)
    receptors <- unique(ctints$interactions)
    nr_receptors <- length(receptors)
    
    return_df <- data.frame(
      celltypes = rep(ctints$celltype[1], length(receptors)),
      annotations = rep("receiver", length(receptors)),
      sigmols = receptors
    )
  }
  
  return(return_df)
}

# for each ct = element in ctints_list[[2]], the expressed l or r are obtained
# returns a list with identical architecture to ctints_list[[2]] 

#-------------------------------------------------------------------------------


extract_lrs_nrs <- function(cpi, ctlrs){
  # cell type pair interactome object
  # ctlrs = info ligands or receptors per CELL TYPE from extract_lrs_info()
  
  print_list <- list()
  
  # for loop because two lists are required (cpi and ctlrs)
  for(j in 1:length(cpi)){
    
    print(j)
    
    emitter_ct <- unique(cpi[[j]]$Celltypes$emitter)
    receiver_ct <- unique(cpi[[j]]$Celltypes$receiver)
    cts <- c(emitter_ct, receiver_ct)
    
    print_df <- data.frame(
      celltypes = vector(mode = "character", length = length(cts)),
      nr_lrs = vector(mode = "numeric", length = length(cts)),
      annotation = vector(mode = "character", length = length(cts))
    )
    
    for(i in 1:length(cts)){
      print_df$celltypes[i] <- ctlrs[[j]][[i]]$celltypes[i]
      print_df$nr_lrs[i] <- nrow(ctlrs[[j]][[i]])
      print_df$annotation[i] <- ctlrs[[j]][[i]]$annotation[i]
    }
    
    print_df <- print_df[order(print_df$nr_lrs, decreasing = TRUE),]
    print_list[[j]] <- print_df
  }
  return(print_list)
}

# returns a df per cpi object with the nr of l or r per cell type
# corresponds to "overview" df of extract_ctp_info() and extract_ct_info()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# DIMENSIONALITY REDUCTION FUNCTIONS 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

get_pca <- function(cci, cts = NULL, sample = NULL, variable = NULL){
  # cts = vector of specific cell types of interest
  # sample = name of a qualitative variable (e.g. Age, Species)
  # variable ? = qualitative variables of each sample (e.g. yng/old for Age)
  
  
  #subset and transpose, keep cts of interest
  if(is.null(cts) == TRUE){
    temp_df <- data.frame(t(cci$Score))
  }else{
    for(i in 1:length(cts))
      temp_df <- data.frame(t(cci$Score[,grep(cts[i], colnames(cci$Score))]))
  }
  
  # if no sample is given calculation is quite straight forward
  if(is.null(sample) == TRUE){
    
    print_pca <- FactoMineR::PCA(temp_df, graph = FALSE)
    
    # if sample is give, information needs to be added to the temp_df
  }else if(is.null(sample) == FALSE){
    
    #add an empty vector to contain samples or other variables
    for(j in 1:length(sample)){
      temp_df[,(ncol(temp_df)+1)] <- vector()
      colnames(temp_df)[ncol(temp_df)] <- sample[j]
      temp_df[,ncol(temp_df)] <- variable[[j]]
    }
    
    # indicate which columns of temp_df are just variables, not numbers
    quali_cols <- which(!is.na(match(colnames(temp_df), sample)))
    print_pca <- FactoMineR::PCA(temp_df, graph = FALSE, quali.sup = quali_cols)
    print(quali_cols)
  }
  
  return(print_pca)
}

# calculates PCA coordinates for cpi objects

#-------------------------------------------------------------------------------
