# this script contains small helper functions that can be used to make the SCE
# objects more convenient to work with

#-------------------------------------------------------------------------------

library(scran, quietly = TRUE) 
library(scater, quietly = TRUE) 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

# SCE factorization for nicer plots

# all cell types

factor_reference_cts <- function(sce){
  
  cell_types_baccin_ordered <- c(
    "LMPPs",
    "Ery/Mk prog.",
    "Mk prog.",
    "Ery prog.",
    "Erythroblasts",
    "Gran/Mono prog.",
    "Eo/Baso prog.",
    "Mono prog.",
    "Monocytes",
    "Dendritic cells",
    "Neutro prog.",
    "Neutrophils",
    "large pre-B.",
    "small pre-B.",
    "pro-B",
    "B cell",
    "NK cells",
    "T cells",
    
    "Ng2+ MSCs",
    "Fibro/Chondro p.",
    "Myofibroblasts",
    "Arteriolar fibro.",
    "Endosteal fibro.",
    "Stromal fibro.",
    "Arteriolar ECs",
    "Sinusoidal ECs",
    "Chondrocytes",
    "Osteoblasts",
    "Smooth muscle",
    "Schwann cells",
    "Adipo-CAR",
    "Osteo-CAR",
    
    "unassigned"
  )
  
  cell_types_dahlin_ordered <- c(
    "HSC LT",
    "HSC ST",
    "MPP1",
    "MPP2",
    "LMPP",
    "MEP",
    "Mkp",
    "Ery", 
    "Erythroblast",
    "Mast",
    "Eosinophil",
    "GMP",
    "Monocyte Progenitor",
    "preNeu",
    "proNeu",
    "Neu",
    "Monocyte",
    "pDC",
    "cDC",
    "Infection related myeloid",
    "CLP",
    "preB",
    "proB",
    "ImmatureB", 
    "Plasma",
    "NK",
    "T",
    "unassigned"
  )
  
  cell_types_dolgalev_ordered <- c(
    "MSPC-Osteo",
    "MSPC-Adipo",
    "Osteo",
    "Osteoblasts",
    "Chondrocytes",
    "Myofibroblasts",
    "Fibroblasts",
    "Pericytes",
    "Smooth-muscle",
    "Schwann-cells",
    "EC-Arteriar",
    "EC-Arteriolar",
    "EC-Sinusoidal",
    "C",
    "unassigned"
  )
  
  cell_types_baccin_ordered_clst <- cell_types_baccin_ordered[
    cell_types_baccin_ordered %in% sce$baccin_celltype_scmapclust]
  sce$baccin_celltype_scmapcell[is.na(sce$baccin_celltype_scmapcell)] <- "unassigned"
  cell_types_baccin_ordered_cell <- cell_types_baccin_ordered[
    cell_types_baccin_ordered %in% sce$baccin_celltype_scmapcell]
  
  sce$baccin_celltype_scmapclust <- factor(sce$baccin_celltype_scmapclust,
                                           levels = cell_types_baccin_ordered_clst)
  sce$baccin_celltype_scmapcell <- factor(sce$baccin_celltype_scmapcell,
                                           levels = cell_types_baccin_ordered_cell)
  
  cell_types_dahlin_ordered_clst <- cell_types_dahlin_ordered[
    cell_types_dahlin_ordered %in% sce$dahlin_celltype_scmapclust]
  sce$dahlin_celltype_scmapcell[is.na(sce$dahlin_celltype_scmapcell)] <- "unassigned"
  cell_types_baccin_ordered_cell <- cell_types_dahlin_ordered[
    cell_types_dahlin_ordered %in% sce$dahlin_celltype_scmapcell]

  sce$dahlin_celltype_scmapclust <- factor(sce$dahlin_celltype_scmapclust,
                                           levels = cell_types_dahlin_ordered_clst)
  sce$dahlin_celltype_scmapcell <- factor(sce$dahlin_celltype_scmapcell,
                                          levels = cell_types_baccin_ordered_cell)

  cell_types_dolgalev_ordered_clst <- cell_types_dolgalev_ordered[
    cell_types_dolgalev_ordered %in% sce$dolgalev_celltype_scmapclust]
  sce$dahlin_celltype_scmapcell[is.na(sce$dahlin_celltype_scmapcell)] <- "unassigned"
  cell_types_dolgalev_ordered_cell <- cell_types_dolgalev_ordered[
    cell_types_dolgalev_ordered %in% sce$dahlin_celltype_scmapcell]
  
  sce$dolgalev_celltype_scmapclust <- factor(sce$dolgalev_celltype_scmapclust,
                                             levels = cell_types_dolgalev_ordered_clst)
  sce$dolgalev_celltype_scmapcell <- factor(sce$dolgalev_celltype_scmapcell,
                                            levels = cell_types_dolgalev_ordered_cell)
  return(sce)
}

#-------------------------------------------------------------------------------

# generic dimensionality reduction for quick UMAP visualisation of SCE objects
# only for intermediate dimensionality reduction (e.g. for reports):
# - preprocessing at sample level
# - objects merged per fraction, species, or age
reduce_dims <- function(sce, nr_hvgs){
  
  hvgs <- scran::modelGeneVar(sce)
  hvgs <- scran::getTopHVGs(hvgs, n=nr_hvgs)
  
  sce <- scater::runPCA(sce, ncomponents=25, subset_row = hvgs) 
  sce <- scater::runUMAP(sce, dimred = 'PCA',
                         external_neighbors=TRUE, subset_row = hvgs)
  return(sce)  
}

#-------------------------------------------------------------------------------
