#-------------------------------------------------------------------------------
# Functions for plotting
#-------------------------------------------------------------------------------

library(cowplot, quietly = TRUE) 
library(ggpubr, quietly = TRUE) 
library(tidyverse, quietly = TRUE) 

#-------------------------------------------------------------------------------
# basic UMAP, must contain a Dimred called "UMAP", must specify factor for color
umap_base <- function(sce, color_by){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["UMAP"]]),
    aes(x = X1,
        y = X2, 
        color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    ggplot2::geom_point(size = 0.01)+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.position = "none", 
                   axis.text = element_blank(),
                   axis.ticks = element_blank())+
    ggplot2::ylab("UMAP 2")+
    ggplot2::xlab("UMAP 1")
  
  return(base)
}

umap_base_l <- function(sce, color_by){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["UMAP"]]),
    aes(x = X1,
        y = X2, 
        color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    ggplot2::geom_point(size = 0.01)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text = element_blank(),
                   axis.ticks = element_blank())+
    ggplot2::ylab("UMAP 2")+
    ggplot2::xlab("UMAP 1")+
    ggplot2::scale_color_discrete(name = color_by)
  
  return(base)
}

#-------------------------------------------------------------------------------
# UMAP for gene expression
umap_gene <- function(sce, color_by){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["UMAP"]]),
    aes(x = X1, 
        y = X2, 
        color = SummarizedExperiment::assays(sce)[["logcounts"]][
          rownames(SummarizedExperiment::assays(sce)[["logcounts"]]) == color_by,]))+
    ggplot2::scale_color_gradientn(color_by, colors = c("black", 
                                                       "darkorange3", 
                                                       "orange", 
                                                       "lightgoldenrod"))+
    ggplot2::geom_point(size = 0.01)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text = element_blank(),
                  axis.ticks = element_blank())+
    ggplot2::ylab("UMAP 2")+
    ggplot2::xlab("UMAP 1")

  return(base)
}

pca_gene <- function(sce, color_by){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["PCA"]]),
    aes(x = PC1, 
        y = PC2, 
        color = SummarizedExperiment::assays(sce)[["logcounts"]][
          rownames(SummarizedExperiment::assays(sce)[["logcounts"]]) == color_by,]))+
    ggplot2::scale_color_gradientn(color_by, colors = c("black", 
                                                       "darkorange3", 
                                                       "orange", 
                                                       "lightgoldenrod"))+
    ggplot2::geom_point(size = 0.3)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text = element_blank(),
                  axis.ticks = element_blank())+
    ggplot2::ylab("UMAP 2")+
    ggplot2::xlab("UMAP 1")
  
  return(base)
}

#-------------------------------------------------------------------------------
# UMAP for (gene expression) programs
# not sure if still in use
umap_program <- function(sce, species, program){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["UMAP"]]),
    aes(x = X1, 
        y = X2, 
        color = metadata(sce)$usage_list[[species]][,program]))+
    ggplot2::scale_color_gradientn(name = base::paste0(program, " usage"),
                                   colors = c("black", 
                                              "darkorange3", 
                                              "orange", 
                                              "lightgoldenrod"))+
    ggplot2::geom_point(size = 0.01)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text = element_blank(),
                   axis.ticks = element_blank())+
    ggplot2::ylab("UMAP 2")+
    ggplot2::xlab("UMAP 1")
  
  return(base)
}

#-------------------------------------------------------------------------------
# UMAP for continuous colData
umap_cont <- function(sce, color_by){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["UMAP"]]),
    aes(x = X1, 
        y = X2, 
        color = colData(sce)[,which(colnames(colData(sce))== color_by)]))+
    ggplot2::scale_color_gradientn(name = color_by, 
                                   colors = c("grey80", "blue"))+
    ggplot2::geom_point(size = 0.001)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text = element_blank(),
                   axis.ticks = element_blank())+
    ggplot2::ylab("UMAP 2")+
    ggplot2::xlab("UMAP 1")
  
  return(base)
}

#-------------------------------------------------------------------------------
# basic UMAP only for legend generation (point sizes visible in legends)
umap_legend <- function(sce, color_by){

  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["UMAP"]]),
    aes(x = X1,
        y = X2, 
        color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    ggplot2::geom_point(size = 2)+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.key.size = unit(0.3, "cm"), 
          legend.key.width = unit(0.3, "cm"),
          legend.spacing = unit(0.06, "cm"), 
          legend.text = element_text(size = 7))
  
  return(base)
  
}

#-------------------------------------------------------------------------------
# basic PCA Dims 1 and 3,  must specify Factor for color
pca_base2 <- function(sce, color_by){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["PCA"]]),
    aes(x = PC1,
        y = PC2, 
        color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    ggplot2::geom_point(size = 0.01)+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.position = "none", axis.text = element_blank())+
    ggplot2::ylab("PCA 2")+
    ggplot2::xlab("PCA 1")
  
  return(base)
  
}

pca_base3 <- function(sce, color_by){
  
  base <- ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["PCA"]]),
    aes(x = PC1, 
        y = PC3, 
        color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    ggplot2::geom_point(size = 0.01)+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.position = "none", axis.text = element_blank())+
    ggplot2::ylab("PCA 3")+
    ggplot2::xlab("PCA 1")
  
  return(base)
  
}

#-------------------------------------------------------------------------------
# basic PCA only for legend generation (point sizes visible in legends)
pca_legend <- function(sce, color_by){
  
  base <- ggplot2::ggplot(
    base::data.frame(SingleCellExperiment::reducedDims(sce)[["PCA"]]),
    aes(x = PC1,
        y = PC2, 
        color = colData(sce)[,colnames(colData(sce)) == color_by]))+
    ggplot2::geom_point(size = 2)+
    ggplot2::theme_classic()+
    ggplot2::theme(legend.key.size = unit(0.3, "cm"), 
                   legend.key.width = unit(0.3, "cm"),
                   legend.spacing = unit(0.06, "cm"), 
                   legend.text = element_text(size = 7))
  
  return(base)
  
}