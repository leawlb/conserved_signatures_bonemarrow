#-------------------------------------------------------------------------------

library(scater)
library(scran)
library(bluster)
library(tidyverse)
set.seed(37)

#-------------------------------------------------------------------------------

core_cons_list <- readRDS(snakemake@input[["core_cons_list"]])

sce <- readRDS(file = snakemake@input[["sce_input"]])

sce_baccin <- readRDS(snakemake@input[["sce_baccin_in"]])
sce_dolgalev <- readRDS(snakemake@input[["sce_dolgalev"]])
sce_dahlin <- readRDS(snakemake@input[["sce_dahlin"]])

k_louvain <- snakemake@params[["k_louvain"]]

fraction_curr <- snakemake@wildcards[["fraction"]]

core_cons_df_list <- lapply(as.list(names(core_cons_list)), function(c){
  
  core_cons_list[[c]]$marker_cons_ct$cell_type <- c
  if(fraction_curr == "hsc"){
    core_cons_list[[c]]$marker_cons_ct$assignment <- "receiver"
  }else if(fraction_curr == "str"){
    core_cons_list[[c]]$marker_cons_ct$assignment <- "emitter"
  }
  
  return(core_cons_list[[c]]$marker_cons_ct)
  
})

core_cons_df <- bind_rows(core_cons_df_list)

cc_list <- unique(core_cons_df$gene[which(core_cons_df$cons_emf == "TRUE")])

#-------------------------------------------------------------------------------

# subset to only relevant genes
# subset baccin to cells that are part of the correct fraction
sce_baccin <- sce_baccin[rownames(sce_baccin) %in% cc_list,]

table(sce_baccin$identity_ref)

if(fraction_curr == "str"){
  sce_ref2 <- sce_dolgalev
  sce_ref2 <- sce_ref2[rownames(sce_ref2) %in% cc_list,]
  
  sce_baccin <- sce_baccin[,sce_baccin$identity_ref %in% c(
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
    "Sinusoidal ECs",
    "Smooth muscle",
    "Stromal fibro.")]
  
  # remove cycling cells and cells that are not in our dataset
  sce_ref2 <- sce_ref2[,!sce_ref2$identity_ref %in% c(
    "Schwann-cells",
    "C")]
  
}else if(fraction_curr == "hsc"){
  sce_ref2 <- sce_dahlin
  sce_ref2 <- sce_ref2[rownames(sce_ref2) %in% cc_list,]
  
  sce_baccin <- sce_baccin[,sce_baccin$identity_ref %in% c(
    "B cell",
    "Dendritic cells",
    "Eo/Baso prog.",
    "Ery prog.",
    "Ery/Mk prog.",
    "Erythroblast",
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
    "small pre-B.")]
  
  # remove really mature cells that are not in our dataset at all
  sce_ref2 <- sce_ref2[,!sce_ref2$identity_ref %in% c(
    
    "Monocyte",
    "preNeu",
    "Neu",
    "cDC",
    "pDC",
    "T",
    "NK",
    "Plasma",
    "Infection related myeloid",
    "Ery",
    "unassigned",
    "Erythroblast",
    "Eosinophil")]
  
}

print(length(cc_list))
print(nrow(sce_baccin))
print(nrow(sce_ref2))

table(sce_baccin$identity_ref)
table(sce_ref2$identity_ref)

#-------------------------------------------------------------------------------

# re-calculate PCA (Without batch correction)
sce_baccin <- runPCA(sce_baccin, ncomponents=25) 
sce_ref2 <- runPCA(sce_ref2, ncomponents=25) 

# calculate clusters
louvain_baccin <- clusterCells(sce_baccin, use.dimred="PCA", assay.type = NULL,
                               BLUSPARAM=NNGraphParam(cluster.fun="louvain", 
                                                      k = k_louvain))

louvain_ref2 <- clusterCells(sce_ref2, use.dimred="PCA", assay.type = NULL,
                             BLUSPARAM=NNGraphParam(cluster.fun="louvain", 
                                                    k = k_louvain))

sce_baccin$cluster_core <- louvain_baccin
sce_ref2$cluster_core <- louvain_ref2

sce_baccin$cluster_core <- factor(sce_baccin$cluster_core,
                                  levels = sort(
                                    unique(sce_baccin$cluster_core))) 
sce_ref2$cluster_core <- factor(sce_ref2$cluster_core,
                                levels = sort(unique(sce_ref2$cluster_core))) 

saveRDS(sce_baccin, file = snakemake@output[["sce_baccin_out"]])
saveRDS(sce_ref2, file = snakemake@output[["sce_ref2_out"]])