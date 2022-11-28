# this script contains small helper functions that can be used to make the SCE
# objects more convenient to work with

#-------------------------------------------------------------------------------

library(scran, quietly = TRUE) 
library(scater, quietly = TRUE) 

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

# generic dimensionality reduction for quick UMAP visualisation of SCE objects
reduce_dims <- function(sce, nr_hvgs){
  
  hvgs <- modelGeneVar(sce)
  hvgs <- getTopHVGs(hvgs, n=nr_hvgs)
  
  sce <- runPCA(sce, ncomponents=25, subset_row = hvgs) 
  sce <- runUMAP(sce, dimred = 'PCA',
                 external_neighbors=TRUE, subset_row = hvgs)
  return(sce)  
}

#-------------------------------------------------------------------------------

