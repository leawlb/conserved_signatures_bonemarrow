# this script contains small helper functions that can be used to make the SCE
# objects more convenient to work with

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Find contamination from Haas reference

hsc_cts_all <- c(
  "B cell",
  "Dendritic cells",
  "Eo/Baso prog.",
  "Ery prog.",
  "Ery/Mk prog.",
  "Erythroblasts",
  "Gran/Mono prog.",
  "LMPPs",
  "Mk prog.",
  "Mono prog.",
  "Monocytes",
  "NK cells",
  "Neutro prog.",
  "Neutrophils",
  "T cells",
  "large pre-B.",
  "pro-B",
  "small pre-B."
)

str_cts_all <- c(
  "Adipo-CAR",
  "Arteriolar ECs",
  "Arteriolar fibro.",
  "Chondrocytes",
  "Endosteal fibro.",
  "Fibro/Chondro p.",
  "Myofibroblasts",
  "Ng2+ MSCs",
  "Osteo-CAR",
  "Osteoblasts",
  "Schwann cells",
  "Sinusoidal ECs",
  "Smooth muscle",
  "Stromal fibro."
)

find_contamination_ref_all <- function(sce, input){
  
  sce$right_fraction <- rep("TRUE", n = ncol(sce))
  
  if(input == "stromal"){
    wrong_pos <- which(sce$Identity_ref_all %in% hsc_cts_all)
    sce$right_fraction[wrong_pos] <- "FALSE"
  
  }else if(input == "hsc"){
    wrong_pos <- which(sce$Identity_ref_all %in% str_cts_all)
    sce$right_fraction[wrong_pos] <- "FALSE"
  }
  return(sce)
}

#-------------------------------------------------------------------------------

# Get Mature cells from Haas reference

hsc_mature_cts_all <- c(
  "B cell",
  "Dendritic cells",
  "Monocytes",
  "small pre-B.",
  "T cells",
  "NK cells",
  "Neutrophils",
  "large pre-B."
)

get_mature_cts_all <- function(sce){
  sce$mature_cells <- rep("FALSE", n = ncol(sce))
  sce$mature_cells[sce$Identity_ref_all %in% hsc_mature_cts_all] <- "TRUE"
  return(sce)
}

#-------------------------------------------------------------------------------

# a vector of ct names in the appropriate order
# use the exact same ct names as intended after changing them
# names should not be able to grep each other or contain + or ()

cts_ordered <- c(
 
  "Ng2-MSCs",#
  "Fibro/Chondro progs",#
  "Myofibroblasts",#
  "Arteriolar fibros",#
  "Endosteal fibros",#
  "Stromal fibros",#
  
  "Adipo-CARs",#
  "Osteo-CARs",#
  
  "Arteriolar ECs",#
  "Sinusoidal ECs",#
  "Smooth muscle",#
  
  "Chondrocytes",#
  "Osteoblasts",#
  
  "Schwann cells",#
  
  "LT-HSCs",#
  "ST-HSCs/MPPs",#
  
  "LMPPs",#
  
  "Small pre-BCs",#
  "Large pre-BCs",#
  "Pro-BCs",#
  "Immature B cells",#
  "Mature B cells",#
  "T cells/NK cells",#
  
  "Common myeloid progs",#
  "Monocyte-primed progs",  #
  "Dendritic-primed progs",  #
  "Neutrophil progs",#
  "Baso/Eo/Mast progs",#
  
  "Monocytes",#
  "Dendritic cells",#
  "Neutrophils",#
  
  "Ery/Mk progs",#
  "Megakaryocyte progs",#
  "Erythrocyte progs",#
  "Erythroids"#
  
)

#-------------------------------------------------------------------------------

# a function that changes the original ct names into the required ct names
# ct names should be stored in a ColData slot called "Identity"
# original ct names will be replaced by required ct names

change_cts <- function(sce){
  
  # if the Identity is already factored, unfactor it
  if(class(colData(sce)$Identity) == "factor"){
    colData(sce)$Identity <- unfactor(colData(sce)$Identity)
  }
  
  # these need to be changed/shortened/summarized
  sce$Identity[colData(sce)$Identity == "Multipotent progenitor (cycling)"] <- "ST-HSCs/MPPs"
  sce$Identity[colData(sce)$Identity == "ST-HSCs"] <- "ST-HSCs/MPPs"

  sce$Identity[colData(sce)$Identity == "Lymphoid primed multipotent progenitor (LMPP)"] <- "LMPPs"
  sce$Identity[colData(sce)$Identity == "Small pre B-cell"] <- "Small pre-BCs"
  sce$Identity[colData(sce)$Identity == "small pre-B."] <- "Small pre-BCs"
  sce$Identity[colData(sce)$Identity == "Large pre B-cell"] <- "Large pre-BCs"
  sce$Identity[colData(sce)$Identity == "large pre-B."] <- "Large pre-BCs"
  sce$Identity[colData(sce)$Identity == "Pro B-cell"] <- "Pro-BCs"
  sce$Identity[colData(sce)$Identity == "pro-B"] <- "Pro-BCs"
  sce$Identity[colData(sce)$Identity == "Immature B-cell"] <- "Immature B cells"
  sce$Identity[colData(sce)$Identity == "Mature B-cell"] <- "Mature B cells"
  sce$Identity[colData(sce)$Identity == "B cell"] <- "Mature B cells"
  sce$Identity[colData(sce)$Identity == "T-cells / NK Cells"] <- "T cells/NK cells"
  sce$Identity[colData(sce)$Identity == "NK cells"] <- "T cells/NK cells"
  sce$Identity[colData(sce)$Identity == "T cells"] <- "T cells/NK cells"
  
  sce$Identity[colData(sce)$Identity == "Common myeloid progenitor"] <- "Common myeloid progs"
  sce$Identity[colData(sce)$Identity == "Basophil / eosinophil / mast cell progenitor ( + mature basophils)"] <- "Baso/Eo/Mast progs"
  sce$Identity[colData(sce)$Identity == "Eo/Baso prog."] <- "Baso/Eo/Mast progs"
  sce$Identity[colData(sce)$Identity == "Neutrophil progenitor"] <- "Neutrophil progs"
  sce$Identity[colData(sce)$Identity == "Neutro prog."] <- "Neutrophil progs"
  
  sce$Identity[colData(sce)$Identity == "Monocyte / Dendritic cell progenitor #1 (more DCs)"] <- "Dendritic-primed progs"
  sce$Identity[colData(sce)$Identity == "Monocyte / Dendritic cell progenitor #2 (more Mono)"] <- "Monocyte-primed progs"
  sce$Identity[colData(sce)$Identity == "Mono prog."] <- "Monocyte-primed progs"
  sce$Identity[colData(sce)$Identity == "Monocyte prog."] <- "Monocyte-primed progs"
  sce$Identity[colData(sce)$Identity == "Gran/Mono prog."] <- "Monocyte-primed progs"
  
  sce$Identity[colData(sce)$Identity == "Neutrophil"] <- "Neutrophils"
  
  sce$Identity[colData(sce)$Identity == "Megakaryocyte / erythroid primed MPP"] <- "Ery/Mk progs"
  sce$Identity[colData(sce)$Identity == "Megakaryocyte / erythroid progenitor"] <- "Ery/Mk progs"
  sce$Identity[colData(sce)$Identity == "Megakaryocyte progenitor"] <- "Megakaryocyte progs"
  sce$Identity[colData(sce)$Identity == "Mk prog."] <- "Megakaryocyte progs"
  sce$Identity[colData(sce)$Identity == "Erythroid Progenitor"] <- "Erythrocyte progs"
  sce$Identity[colData(sce)$Identity == "Ery prog."] <- "Erythrocyte progs"
  sce$Identity[colData(sce)$Identity == "Erythroblasts"] <- "Erythroids"
  sce$Identity[colData(sce)$Identity == "Erythroid"] <- "Erythroids"
  sce$Identity[colData(sce)$Identity == "Ery/Mk prog."] <- "Ery/Mk progs"

  sce$Identity[colData(sce)$Identity == "Ng2+ MSCs"] <- "Ng2-MSCs"
  sce$Identity[colData(sce)$Identity == "Fibro/Chondro p."] <- "Fibro/Chondro progs"
  sce$Identity[colData(sce)$Identity == "Adipo-CAR"] <- "Adipo-CARs"
  sce$Identity[colData(sce)$Identity == "Osteo-CAR"] <- "Osteo-CARs"
  sce$Identity[colData(sce)$Identity == "Arteriolar fibro."] <- "Arteriolar fibros"
  sce$Identity[colData(sce)$Identity == "Endosteal fibro."] <- "Endosteal fibros"
  sce$Identity[colData(sce)$Identity == "Stromal fibro."] <- "Stromal fibros"
  
  # remove all unclear/dying cells
  if("Unclear / dying cells" %in% colData(sce)$Identity){
    sce <-  sce[,-which(colData(sce)$Identity == "Unclear / dying cells")]
  }
  # factorize 
  colData(sce)$Identity <- factor(
    colData(sce)$Identity, levels = cts_ordered[
      which(cts_ordered %in% unique(colData(sce)$Identity))]
  )
  return(sce)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# a function to add shared ENSMUS IDs from insp_ref in correct positions
# requires loading of insp_ref_shared_man !!!!

get_shared_ensids <- function(sce){
  
  # extract the species to know which IDs to match
  species <- colData(sce)$Species[1]
  # mmus is called "mbl6" in insp_ref
  if(species == "mmus"){
    species <- "mbl6"
  }
  
  # find all intersecting genes (since there may be duplicates in insp_ref)
  ints_geneids <- intersect(insp_ref$id[,grep(species, colnames(insp_ref$id))], 
                            rowData(sce)$ID)
  
  # mark the positions in sce that fit the ints_geneids = shared with insp_ref
  match_pos_sce <- match(ints_geneids, rowData(sce)$ID)
  
  # make a slot to store whether transcripts are shared with all other species
  rowData(sce)$Shared <- vector(length = nrow(sce))
  rowData(sce)$Shared[match_pos_sce] <- "shared"
  
  # get the positions in insp_ref that fit the positions in sce
  match_pos_insp <- match(rowData(sce)$ID[rowData(sce)$Shared == "shared"],
                          insp_ref$id[,grep(species, colnames(insp_ref$id))])
  
  # make new slot for ENSMUS-IDs and fill
  rowData(sce)$ENSMUS_ID <- vector(length = nrow(sce))
  rowData(sce)$ENSMUS_ID[
    rowData(sce)$Shared == "shared"] <- insp_ref$id$ensembl_mbl6_id[match_pos_insp]
  
  return(sce)
}

