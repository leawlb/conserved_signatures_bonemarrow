#-------------------------------------------------------------------------------
# Functions to make my life regarding plots easier
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
# basic UMAP, must contain a Dimred called "UMAP", must specify Factor for color
umap_base <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["UMAP"]]),
                 aes(x = X1, y = X2, color = color_by))+
    geom_point(size = 0.01)+
    theme_classic()+
    theme(legend.position = "none", axis.text = element_blank())+
    ylab("UMAP 2")+
    xlab("UMAP 1")
  
  return(base)
  
}

#-------------------------------------------------------------------------------
# basic UMAP only for legend generation (point sizes visible in legends)
umap_legend <- function(sce, color_by){
  
  base <- ggplot(data.frame(reducedDims(sce)[["UMAP"]]),
                 aes(x = X1, y = X2, color = sce$Identity_ref_all))+
    geom_point(size = 2)+
    theme_classic+
    theme(legend.key.size = unit(0.3, "cm"), 
          legend.key.width = unit(0.3, "cm"),
          legend.spacing = unit(0.06, "cm"), 
          legend.text = element_text(size = 7))
  
  return(base)
  
}