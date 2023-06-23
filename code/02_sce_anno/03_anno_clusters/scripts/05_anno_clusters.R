#-------------------------------------------------------------------------------
library(DropletUtils, quietly = TRUE) 

#TODO: either completely delete when not needed anymore or take from CSV

sce <- readRDS(file = snakemake@input[["sce_input"]])
fraction_curr <- sce$Fraction_ID[1]

# manually assign first layer of cell type annotations

sce$prelimanno_broad <- vector(length = ncol(sce))

if(fraction_curr == "hsc"){
  
  sce$prelimanno_broad[sce$cluster_louvain == "1"] <- "Myeloid/Granu early"
  sce$prelimanno_broad[sce$cluster_louvain == "2"] <- "MPPs early"
  sce$prelimanno_broad[sce$cluster_louvain == "3"] <- "Erythroid"
  sce$prelimanno_broad[sce$cluster_louvain == "4"] <- "HSCs"
  sce$prelimanno_broad[sce$cluster_louvain == "5"] <- "MPPs differentiating"
  sce$prelimanno_broad[sce$cluster_louvain == "6"] <- "Myeloid/Granu late"
  sce$prelimanno_broad[sce$cluster_louvain == "7"] <- "Cycling1"
  sce$prelimanno_broad[sce$cluster_louvain == "8"] <- "Baso/Eo/Mast progs."
  sce$prelimanno_broad[sce$cluster_louvain == "9"] <- "Mk/Erythroid progs."
  sce$prelimanno_broad[sce$cluster_louvain == "10"] <- "Mature cells, other"
  sce$prelimanno_broad[sce$cluster_louvain == "11"] <- "Mk"
  sce$prelimanno_broad[sce$cluster_louvain == "12"] <- "Cycling2"
  sce$prelimanno_broad[sce$cluster_louvain == "13"] <- "Mature cells, myeloid"
  
  sce$prelimanno_broad <- factor(sce$prelimanno_broad,
                                     levels = c(
                                       "HSCs",
                                       "MPPs early",
                                       "MPPs differentiating",
                                       "Myeloid/Granu early",
                                       "Myeloid/Granu late",
                                       "Baso/Eo/Mast progs.",
                                       "Mk",
                                       "Mk/Erythroid progs.",
                                       "Erythroid",
                                       "Cycling1",
                                       "Cycling2",
                                       "Mature cells, myeloid",
                                       "Mature cells, other"
                                     ))

}else if(fraction_curr == "str"){
  
  sce$prelimanno_broad[sce$cluster_louvain == "1"] <- "ery lineage"
  sce$prelimanno_broad[sce$cluster_louvain == "2"] <- "mono lineage 1"
  sce$prelimanno_broad[sce$cluster_louvain == "3"] <- "mast/eo or granu"
  sce$prelimanno_broad[sce$cluster_louvain == "4"] <- "Mesenchymal"
  sce$prelimanno_broad[sce$cluster_louvain == "5"] <- "mono lineage 2"
  sce$prelimanno_broad[sce$cluster_louvain == "6"] <- "Arteriolar/Mixed ECs"
  sce$prelimanno_broad[sce$cluster_louvain == "7"] <- "neutro lineage"
  sce$prelimanno_broad[sce$cluster_louvain == "8"] <- "Fibroblasts"
  sce$prelimanno_broad[sce$cluster_louvain == "9"] <- "CAR cells 1"
  sce$prelimanno_broad[sce$cluster_louvain == "10"] <- "Arteriolar/Capillary ECs"
  sce$prelimanno_broad[sce$cluster_louvain == "11"] <- "lymphoid 1"
  sce$prelimanno_broad[sce$cluster_louvain == "12"] <- "Pericytes/Smooth muscle cells"
  sce$prelimanno_broad[sce$cluster_louvain == "13"] <- "Sinusoidal ECs"
  sce$prelimanno_broad[sce$cluster_louvain == "14"] <- "CAR cells 2"
  sce$prelimanno_broad[sce$cluster_louvain == "15"] <- "lymphoid 2"
  sce$prelimanno_broad[sce$cluster_louvain == "16"] <- "Skeletal muscle cells"
  sce$prelimanno_broad[sce$cluster_louvain == "17"] <- "lymphoid 3"
  
  sce$prelimanno_broad <- factor(sce$prelimanno_broad,
                                 levels = c(
                                   "Mesenchymal",
                                   "Fibroblasts",
                                   "CAR cells 1",
                                   "CAR cells 2",
                                   "Arteriolar/Capillary ECs",
                                   "Arteriolar/Mixed ECs",
                                   "Sinusoidal ECs",
                                   "Pericytes/Smooth muscle cells",
                                   "Skeletal muscle cells",
                                   "mono lineage 1",
                                   "mono lineage 2",
                                   "neutro lineage",
                                   "lymphoid 1",
                                   "lymphoid 2",
                                   "lymphoid 3",
                                   "ery lineage",
                                   "mast/eo or granu"
                                 ))
}

saveRDS(sce, file = snakemake@output[["sce_output"]])
