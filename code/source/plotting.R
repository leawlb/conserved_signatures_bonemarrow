#-------------------------------------------------------------------------------
# Functions to make my life regarding plots easier
#-------------------------------------------------------------------------------

library(cowplot, quietly = TRUE) 
library(ggpubr, quietly = TRUE) 

#-------------------------------------------------------------------------------
# basic UMAP, must contain a Dimred called "UMAP", must specify Factor for color
umap_base <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["UMAP"]]),
                 aes(x = X1, y = X2, 
                     color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    geom_point(size = 0.01)+
    theme_classic()+
    theme(legend.position = "none", axis.text = element_blank())+
    ylab("UMAP 2")+
    xlab("UMAP 1")
  
  return(base)
  
}

umap_base_l <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["UMAP"]]),
                 aes(x = X1, y = X2, 
                     color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    geom_point(size = 0.01)+
    theme_classic()+
    theme(axis.text = element_blank())+
    ylab("UMAP 2")+
    xlab("UMAP 1")+
    scale_color_discrete(name = color_by)
  
  return(base)
  
}

#-------------------------------------------------------------------------------
# UMAP for gene expression
umap_gene <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["UMAP"]]),
                 aes(x = X1, y = X2, 
                     color = assays(sce)[["logcounts"]][rownames(assays(sce)[["logcounts"]]) == color_by,]))+
    scale_color_gradientn(color_by, colors = c("black", "darkorange3", "orange", "lightgoldenrod"))+
    geom_point(size = 0.01)+
    theme_classic()+
    theme(axis.text = element_blank())+
    ylab("UMAP 2")+
    xlab("UMAP 1")

  return(base)
}

#-------------------------------------------------------------------------------
# UMAP for gene expression
umap_program <- function(sce, species, program){
  
  base <- ggplot(data.frame(reducedDims(sce)[["UMAP"]]),
                 aes(x = X1, y = X2, 
                     color = metadata(sce)$usage_list[[species]][,program]))+
    scale_color_gradientn(name = paste0(program, " usage"), 
                          colors = c("black", "darkorange3", 
                                     "orange", "lightgoldenrod"))+
    geom_point(size = 0.01)+
    theme_classic()+
    theme(axis.text = element_blank())+
    ylab("UMAP 2")+
    xlab("UMAP 1")
  
  return(base)
}

#-------------------------------------------------------------------------------
# basic UMAP only for legend generation (point sizes visible in legends)
umap_legend <- function(sce, color_by){

  base <- ggplot(data.frame(reducedDims(sce)[["UMAP"]]),
                 aes(x = X1, y = X2, 
                     color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    geom_point(size = 2)+
    theme_classic()+
    theme(legend.key.size = unit(0.3, "cm"), 
          legend.key.width = unit(0.3, "cm"),
          legend.spacing = unit(0.06, "cm"), 
          legend.text = element_text(size = 7))
  
  return(base)
  
}


#-------------------------------------------------------------------------------
# basic PCA Dims 1 and 3,  must specify Factor for color
pca_base2 <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["PCA"]]),
                 aes(x = PC1, y = PC2, 
                     color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    geom_point(size = 0.01)+
    theme_classic()+
    theme(legend.position = "none", axis.text = element_blank())+
    ylab("PCA 2")+
    xlab("PCA 1")
  
  return(base)
  
}

pca_base3 <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["PCA"]]),
                 aes(x = PC1, y = PC3, 
                     color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    geom_point(size = 0.01)+
    theme_classic()+
    theme(legend.position = "none", axis.text = element_blank())+
    ylab("PCA 3")+
    xlab("PCA 1")
  
  return(base)
  
}

#-------------------------------------------------------------------------------
# basic PCA only for legend generation (point sizes visible in legends)
pca_legend <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["PCA"]]),
                 aes(x = PC1, y = PC2, 
                     color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    geom_point(size = 2)+
    theme_classic()+
    theme(legend.key.size = unit(0.3, "cm"), 
          legend.key.width = unit(0.3, "cm"),
          legend.spacing = unit(0.06, "cm"), 
          legend.text = element_text(size = 7))
  
  return(base)
  
}