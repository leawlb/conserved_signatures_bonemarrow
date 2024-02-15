library(tidyverse)

#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# functions for extraction relevant info from cell type interactomes


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------

## MAIN FUNCTION: STEP 1
## extract important metrics on interactions at IDENTITIY PAIR level:
## EXTRACT IDENTITY PAIR INFO

# returns a list with an overview df on all identity pairs (nr cells, nr ints)
# and a list of large dataframes for each idp (nr cells, nr ints, ranks)

extract_ident_pair_info <- function(cci, sce){
  
  # purpose: extract info on nr of cells, nr of interactions per identity pair
  
  cci <- cci      # CCI object from make_interaction_list()
  sce <- sce      # SCE object
  
  # for each identity pair (idp) i in the cci object 
  # make a df containing info on all interactions
  return_list <- lapply(as.list(1:ncol(cci$Score)), sce = sce, cci = cci, 
                        function(i, sce, cci){
                          
    cci <- cci
    sce <- sce
    
    temp_df <- data.frame(
      "ident_pair" = rep(colnames(cci$Score)[i], nrow(cci$Score)),
      "interaction" = rownames(cci$Score),
      "score" = cci$Score[,i],
      "emitter" = rep(cci$Identities$emitter[i], nrow(cci$Score)),
      "nr_cells_emitter" = vector(length = nrow(cci$Score)),
      "ligandrank" = cci$Ligandrank[,i],
      "receiver" = rep(cci$Identities$receiver[i], nrow(cci$Score)),
      "nr_cells_receiver" = vector(length = nrow(cci$Score)),
      "receptorrank" = cci$Receptorrank[,i],
      "species" = rep(sce$Species_ID[1], nrow(cci$Score)),
      "age" = rep(sce$Age_ID[1], nrow(cci$Score)),
      "condition" = vector(length = nrow(cci$Score))
    )
    
    temp_df$nr_cells_emitter <- length(grep(temp_df$emitter[1], sce$Identity))
    temp_df$nr_cells_receiver <- length(grep(temp_df$receiver[1], sce$Identity))
    
    if(!"condition" %in% colnames(cci$Identities)){
      temp_df$condition <- "main"
    }else if("condition" %in% colnames(cci$Identities)){
      temp_df$condition <- cci$Identities$condition[1]
    }
    
    # remove all non-detected interactions, then order from high to low score
    temp_df <- temp_df[!is.na(temp_df$score),]
    temp_df <- temp_df[order(temp_df$score, decreasing = TRUE),]
    return(temp_df)
    
  })
  names(return_list) <- colnames(cci$Score)
  
  # from each temp_df of each idp i
  # add the most important info into a smaller df for later overview of all idps
  return_df_list <- lapply(return_list, cci = cci, function(temp_df, cci){
    
    temp_df <- temp_df
    cci <- cci
    
    # make empty df
    return_df <- data.frame(
      "ident_pairs" = temp_df$ident_pair[1]
    )

    # enter one row per idp into a large df collecting the nr of ints per idp
    # for better overview
    return_df$nr_ints[1] <- nrow(temp_df) 
    return_df$emitter[1] <- temp_df$emitter[1]
    return_df$nr_cells_emitter[1] <- temp_df$nr_cells_emitter[1]
    return_df$receiver[1] <- temp_df$receiver[1]
    return_df$nr_cells_receiver[1] <- temp_df$nr_cells_receiver[1]
    return_df$species[1] <- temp_df$species[1]
    return_df$age[1] <- temp_df$age[1]
    return_df$condition[1] <- temp_df$condition[1]
    
    return(return_df)
  })
  
  # bind one-row dfs from df list into one overview df on all idps
  names(return_df_list) <- colnames(cci$Score)
  return_df <- bind_rows(return_df_list)
  return_df <- return_df[order(return_df$nr_ints, decreasing = TRUE),]

  ipi_list <- list("overview" = return_df, "separate" = return_list)
  return(ipi_list)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## MAIN FUNCTION: STEP 2
## extract metrics for each IDENTITY (nr interactions)
## EXTRACT IDENTITY INFO

# returns a list with an overview df on all identities (nr cells, nr ints)
# and a list of large dataframes for each identity (nr cells, nr ints, ranks)

extract_ident_info <- function(ipi_list){
  
  # purpose: extract info on identities from identity pair info
  
  ipi_list <- ipi_list      # info on idps from extract_ident_pair_info()
  
  # put all info from ipi_lists into one big dataframe
  ipi_total <- bind_rows(ipi_list[[2]])

  # rearrange list by identities, first emitters
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
  # ipi_total_list now contains ident_pair_infos sorted by emitter/receiver
  # instead of identity pair
  
  # now use this data to make an idi_list with all relevant information
  ident_info_list <- lapply(ipi_total_list, function(ipi_total){
    
    ident_info <- data.frame(
      interaction = unique(ipi_total$interaction),
      identity = rep(ipi_total$identity[1], 
                     length = length(unique(ipi_total$interaction))),
      assignment = rep(ipi_total$assignment[1],
                       length = length(unique(ipi_total$interaction))),
      species = rep(ipi_total$species[1],
                       length = length(unique(ipi_total$interaction))),
      age = rep(ipi_total$age[1],
                       length = length(unique(ipi_total$interaction))),
      condition = rep(ipi_total$condition[1],
                       length = length(unique(ipi_total$interaction)))
      
    )
    
    # add the relevant rank
    if(ident_info$assignment[1] == "emitter"){
      ident_info$rank <-  ipi_total$ligandrank[
        match(ident_info$interaction, ipi_total$interaction)]
      ident_info$nr_cells <- ipi_total$nr_cells_emitter[
        match(ident_info$interaction, ipi_total$interaction)]
    }else if(ident_info$assignment[1] == "receiver"){
      ident_info$rank <-  ipi_total$receptorrank[
        match(ident_info$interaction, ipi_total$interaction)]
      ident_info$nr_cells <- ipi_total$nr_cells_receiver[
        match(ident_info$interaction, ipi_total$interaction)]
    }
    
    ident_info <- separate(ident_info, col = "interaction", 
                           into = c("ligand", "receptor"), sep = "_",
                           remove = FALSE)
    
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
    overview_df$nr_cells[i] <- ident_info_list[[i]]$nr_cells[1]
    overview_df$species[i] <- ident_info_list[[i]]$species[1]
    overview_df$age[i] <- ident_info_list[[i]]$age[1]
    overview_df$condition[i] <- ident_info_list[[i]]$condition[1]
  }
  overview_df <- overview_df[order(overview_df$nr_ints, decreasing = TRUE),]
  print(overview_df)
  
  idi_list <- list(overview = overview_df, separate = ident_info_list)
  return(idi_list)
}

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## MAIN FUNCTION: STEP 3
## extract metrics on ligands or receptors per IDENTITY:
## EXTRACT LIGAND OR RECEPTOR INFO

# returns for one item of an idi_list a new df with info on signaling molecules
# = ligands or receptors, depending on assignment

extract_lrs_info <- function(idi_list){
  
  idi_list <- idi_list # info on idents from extract_ident_info()
  
  if(idi_list$assignment[1] == "emitter"){
    
    ligands <- unique(idi_list$ligand)
    
    return_df <- data.frame(
      identity = rep(idi_list$identity[1], length(ligands)),
      assignment = rep("emitter", length(ligands)),
      sigmol = ligands,
      rank = vector(length = length(ligands)),
      nr_cells = rep(idi_list$nr_cells[1], length(ligands)),
      species = rep(idi_list$species[1], length(ligands)),
      age = rep(idi_list$age[1], length(ligands)),
      condition = rep(idi_list$condition[1], length(ligands))
    )
    
    return_df$rank <- idi_list$rank[match(return_df$sigmol, idi_list$ligand)]
    
  }else if(idi_list$assignment[1] == "receiver"){
    receptors <- unique(idi_list$receptor)
    
    return_df <- data.frame(
      identity = rep(idi_list$identity[1], length(receptors)),
      assignment = rep("receiver", length(receptors)),
      sigmol = receptors,
      rank = vector(length = length(receptors)),
      nr_cells = rep(idi_list$nr_cells[1], length(receptors)),
      species = rep(idi_list$species[1], length(receptors)),
      age = rep(idi_list$age[1], length(receptors)),
      condition = rep(idi_list$condition[1], length(receptors))
    )
    
    return_df$rank <- idi_list$rank[match(return_df$sigmol, idi_list$receptor)]
  }
  return(return_df)
}

## NOTE not all ranks are included!
## = sigmols that are only part "lonely" interactions are not included!
## These are not part of the interactome of a given identity, even though
## they may be highly expressed in those cells.

#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------

## MAIN FUNCTION STEP 4: 
## extract metrics on nr of signaling molecules/ligands/receptors per IDENTITY: 
## EXTRACT RECEPTOR/LIGAND NUMBERS

# returns a df per object with the nr of l or r per cell type
# corresponds to "overview" dfs from identity and identity pair level lists

extract_lrs_nrs <- function(cci, lrs_list){
  
  cci <- cci                # CCI object
  lrs_list <- lrs_list      # lrs_list = info on l/r from extract_lrs_info()

  emitter_ct <- unique(cci$Identities$emitter)
  receiver_ct <- unique(cci$Identities$receiver)
  idents <- c(emitter_ct, receiver_ct)
  
  return_df <- data.frame(
    identity = vector(mode = "character", length = length(idents)),
    nr_lrs = vector(mode = "numeric", length = length(idents)),
    assignment = vector(mode = "character", length = length(idents))
  )
  
  for(i in 1:length(idents)){
    
    return_df$identity[i] <- lrs_list[[i]]$identity[i]
    return_df$assignment[i] <- lrs_list[[i]]$assignment[i]
    return_df$nr_lrs[i] <- nrow(lrs_list[[i]])
    return_df$nr_cells[i] <- lrs_list[[i]]$nr_cells[i]
    return_df$species[i] <- lrs_list[[i]]$species[i]
    return_df$age[i] <- lrs_list[[i]]$age[i]
    return_df$condition[i] <- lrs_list[[i]]$condition[i]
  }
  
  return_df <- return_df[order(return_df$nr_lrs, decreasing = TRUE),]
  
  return(return_df)
}
