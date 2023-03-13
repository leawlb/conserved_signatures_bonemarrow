#-------------------------------------------------------------------------------
library(DropletUtils, quietly = TRUE) 

#TODO: either completely delete when not needed anymore or take from CSV

sce <- readRDS(file = snakemake@input[["sce_10"]])
fraction_curr <- sce$Fraction_ID[1]

# manually assign first layer of cell type annotations

sce$prelimanno_broad <- vector(length = ncol(sce))

if(fraction_curr == "hsc"){
  
  sce$prelimanno_broad[sce$cluster_louvain == "1"] <- "myeloid prog."
  sce$prelimanno_broad[sce$cluster_louvain == "2"] <- "non-cycling MPPs"
  sce$prelimanno_broad[sce$cluster_louvain == "3"] <- "ery prog./blasts"
  sce$prelimanno_broad[sce$cluster_louvain == "4"] <- "HSCs/MPPs"
  sce$prelimanno_broad[sce$cluster_louvain == "5"] <- "active MPPs"
  sce$prelimanno_broad[sce$cluster_louvain == "6"] <- "committed granu prog."
  sce$prelimanno_broad[sce$cluster_louvain == "7"] <- "cycling1"
  sce$prelimanno_broad[sce$cluster_louvain == "8"] <- "baso/eo/mast prog."
  sce$prelimanno_broad[sce$cluster_louvain == "9"] <- "mk/ery prog."
  sce$prelimanno_broad[sce$cluster_louvain == "10"] <- "mature cells, other"
  sce$prelimanno_broad[sce$cluster_louvain == "11"] <- "mk prog."
  sce$prelimanno_broad[sce$cluster_louvain == "12"] <- "cycling2"
  sce$prelimanno_broad[sce$cluster_louvain == "13"] <- "mature cells, myeloid"
  
  sce$prelimanno_broad <- factor(sce$prelimanno_broad,
                                     levels = c(
                                       "HSCs/MPPs",
                                       "non-cycling MPPs",
                                       "active MPPs",
                                       "myeloid prog.",
                                       "committed granu prog.",
                                       "baso/eo/mast prog.",
                                       "mk prog.",
                                       "mk/ery prog.",
                                       "ery prog./blasts",
                                       "cycling1",
                                       "cycling2",
                                       "mature cells, myeloid",
                                       "mature cells, other"
                                     ))

}else if(fraction_curr == "str"){
  
  sce$prelimanno_broad[sce$cluster_louvain == "1"] <- "ery lineage"
  sce$prelimanno_broad[sce$cluster_louvain == "2"] <- "mono lineage early"
  sce$prelimanno_broad[sce$cluster_louvain == "3"] <- "mast/eo or granu"
  sce$prelimanno_broad[sce$cluster_louvain == "4"] <- "early mesenchymal"
  sce$prelimanno_broad[sce$cluster_louvain == "5"] <- "mono lineage late"
  sce$prelimanno_broad[sce$cluster_louvain == "6"] <- "endos1"
  sce$prelimanno_broad[sce$cluster_louvain == "7"] <- "neutro lineage"
  sce$prelimanno_broad[sce$cluster_louvain == "8"] <- "fibroblasts"
  sce$prelimanno_broad[sce$cluster_louvain == "9"] <- "CAR cells"
  sce$prelimanno_broad[sce$cluster_louvain == "10"] <- "endos2"
  sce$prelimanno_broad[sce$cluster_louvain == "11"] <- "lymphoid1"
  sce$prelimanno_broad[sce$cluster_louvain == "12"] <- "pericyte/smooth muscle?"
  sce$prelimanno_broad[sce$cluster_louvain == "13"] <- "endos3"
  sce$prelimanno_broad[sce$cluster_louvain == "14"] <- "late CAR"
  sce$prelimanno_broad[sce$cluster_louvain == "15"] <- "lymphoid2"
  sce$prelimanno_broad[sce$cluster_louvain == "16"] <- "smooth muscle/pericyte?"
  sce$prelimanno_broad[sce$cluster_louvain == "17"] <- "lymphoid3"
  
  sce$prelimanno_broad <- factor(sce$prelimanno_broad,
                                 levels = c(
                                   "early mesenchymal",
                                   "CAR cells",
                                   "late CAR",
                                   "fibroblasts",
                                   "endos1",
                                   "endos2",
                                   "endos3",
                                   "smooth muscle/pericyte?",
                                   "pericyte/smooth muscle?",
                                   "mono lineage early",
                                   "mono lineage late",
                                   "neutro lineage",
                                   "lymphoid1",
                                   "lymphoid2",
                                   "lymphoid3",
                                   "ery lineage",
                                   "mast/eo or granu"
                                 ))
}

saveRDS(sce, file = snakemake@output[["sce_11"]])


