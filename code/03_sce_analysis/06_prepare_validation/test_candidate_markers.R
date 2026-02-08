
source("repositories/Interspecies_BM_phd/code/source/plotting.R")
source("repositories/Interspecies_BM_phd/code/figures/determine_params.R")
colors_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/colors/colors.txt"
source("repositories/Interspecies_BM_phd/code/source/colors.R")

sce_str <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_str.rds")
sce_str_mmus <- sce_str[,sce_str$Species_ID == "mmus"]


#------------------------------------------------------------------------------
# SEC

list_genes_SEC <- c(

  # to be tested
  "Slco2a1",
  "Il6st",

  # positive pan EC
  "Esam",
  "Eng",
  #"Notch1",

  # positive SEC specific
  "Selp",
  "Vcam1",
  "Icam1",
  "Ctnnal1",
  "Cd55",
  #"Tgfbr2",
  "Bmp4",


  # negative AEC specific
  #"Ly6a",
  #"Cd34",
  "Tcf15",
  #"Kitl"

  # NON-TARGET controls
  "Plvap",
  "Stab2",
  "Lrg1"
)

# umaps

lapply(as.list(list_genes_SEC), function(gene){
  umap_gene(sce_str_mmus, gene)
})

# violin

sce_ec <- sce_str_mmus[rownames(sce_str_mmus) %in% c(
  list_genes_SEC
),]

# make a df for visualisation
logcounts_df_ec <- logcounts(sce_ec) %>%
  as.matrix() %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "genes") %>%
  tidyr::pivot_longer(
    cols = colnames(sce_ec),
    names_to = "cell",
    values_to = "logcounts")

# add cluster info
logcounts_df_ec$celltypes <- sce_ec$celltypes[
  match(
    logcounts_df_ec$cell,
    rownames(colData(sce_ec)))]

logcounts_df_ec$genes <- factor(
  logcounts_df_ec$genes, 
  levels = c(list_genes_SEC))

vioplot_ec <- logcounts_df_ec %>%
  ggplot2::ggplot(
    aes(
      x = celltypes,
      y = logcounts,
      color = celltypes))+
    ggplot2::geom_violin(
      scale = "width",
      draw_quantiles = c(0.25, 0.5, 0.75))+
    #ggplot2::scale_y_continuous(
    #  limits = c(0, 6),
    #  breaks = c(0, 5)
    #)+
    theme_all+
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90, 
        hjust = 0.5,
        vjust = 0,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face),
      strip.text.y.right = element_text(
        angle = 0, 
        hjust = 0,
        vjust = 0.5,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face
      ),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none")+
    ggplot2::facet_grid(rows = vars(genes))+
    ggplot2::scale_color_manual("Cell types", values = col_cts_str)

vioplot_ec

# I would like to go with
# - Esam
# - Eng
# - Selp
# - Vcam1 or Icam1
# - Ctnnal1 observed in SECs
# - Tcf15 if negative marker makes sense

# violin SMALL
list_genes_SEC_small <- c(
  "Il6st",
  "Slco2a1",
  "Esam",
  "Eng",
  "Selp",
  "Vcam1",
  "Ctnnal1", 
  "Tcf15" 
)


# umaps

lapply(as.list(list_genes_SEC), function(gene){
  umap_gene(sce_str_mmus, gene)
})

# small violin

sce_ec <- sce_str_mmus[rownames(sce_str_mmus) %in% c(
  list_genes_SEC_small
),]

# make a df for visualisation
logcounts_df_ec <- logcounts(sce_ec) %>%
  as.matrix() %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "genes") %>%
  tidyr::pivot_longer(
    cols = colnames(sce_ec),
    names_to = "cell",
    values_to = "logcounts")

# add cluster info
logcounts_df_ec$celltypes <- sce_ec$celltypes[
  match(
    logcounts_df_ec$cell,
    rownames(colData(sce_ec)))]

logcounts_df_ec$genes <- factor(
  logcounts_df_ec$genes, 
  levels = c(list_genes_SEC_small))

vioplot_ec <- logcounts_df_ec %>%
  ggplot2::ggplot(
    aes(
      x = celltypes,
      y = logcounts,
      color = celltypes))+
    ggplot2::geom_violin(
      scale = "width",
      draw_quantiles = c(0.25, 0.5, 0.75))+
    #ggplot2::scale_y_continuous(
    #  limits = c(0, 6),
    #  breaks = c(0, 5)
    #)+
    theme_all+
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90, 
        hjust = 0.5,
        vjust = 0,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face),
      strip.text.y.right = element_text(
        angle = 0, 
        hjust = 0,
        vjust = 0.5,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face
      ),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none")+
    ggplot2::facet_grid(rows = vars(genes))+
    ggplot2::scale_color_manual("Cell types", values = col_cts_str)

vioplot_ec

#------------------------------------------------------------------------------
# Fibro (TrkB)

list_genes_fibro <- c(

  # to be tested
  "Ntrk2",

  # positive fibro-pop specific
  #"Cxcl14", # very good overlap but adipo CAr is higher
  "Col5a3",
  "Col14a1",
  #"Angpt1", # very good overlap but adipo CAR is higher
  "Col6a3",
  "Adamts12",
  "Icam1",
  #"Nid1",
  #"Hmcn2",
  #"Fndc1",
  #"Ccn3",
  #"Lepr",
  #"Fbln2",
  #"Fbn1",
  #"Rarres2",
  #"Fn1",
  #"Ecm1",
  #"Emilin2"

  # general fibro
  "Ly6a",
  "Cd34",
  "Pdgfra",
  "Dcn",

  # NON-TARGET controls
  "Abca8a",
  "Abca8b",
  "Celf2",
  "Podn",
  "Cygb",
  "Dpt",
  "Ltbp4",
  "Spry2"
  #"Lama2",
  #"Lamb2",
  #"Abi3bp",
  #"Clec3b"
)

# umaps

lapply(as.list(list_genes_fibro), function(gene){
  umap_gene(sce_str_mmus, gene)
})

# violin

sce_fibro <- sce_str_mmus[rownames(sce_str_mmus) %in% c(
  list_genes_fibro
),]

# make a df for visualisation
logcounts_df_fibro <- logcounts(sce_fibro) %>%
  as.matrix() %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "genes") %>%
  tidyr::pivot_longer(
    cols = colnames(sce_fibro),
    names_to = "cell",
    values_to = "logcounts")

# add cluster info
logcounts_df_fibro$celltypes <- sce_fibro$celltypes[
  match(
    logcounts_df_fibro$cell,
    rownames(colData(sce_fibro)))]

logcounts_df_fibro$genes <- factor(
  logcounts_df_fibro$genes, 
  levels = c(list_genes_fibro))

vioplot_fibro <- logcounts_df_fibro %>%
  ggplot2::ggplot(
    aes(
      x = celltypes,
      y = logcounts,
      color = celltypes))+
    ggplot2::geom_violin(
      scale = "width",
      draw_quantiles = c(0.25, 0.5, 0.75))+
    #ggplot2::scale_y_continuous(
    #  limits = c(0, 6),
    #  breaks = c(0, 5)
    #)+
    theme_all+
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90, 
        hjust = 0.5,
        vjust = 0,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face),
      strip.text.y.right = element_text(
        angle = 0, 
        hjust = 0,
        vjust = 0.5,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face
      ),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none")+
    ggplot2::facet_grid(rows = vars(genes))+
    ggplot2::scale_color_manual("Cell types", values = col_cts_str)

vioplot_fibro

# I would like to go with
# - Ly6a
# - Cd34
# - Pdgfra
# - Col5a3 fibroblast subset
# - Col14a1 fibroblast subset


# violin small

list_genes_fibro_small <- c(
  "Ntrk2",
  "Ly6a",
  "Cd34",
  "Pdgfra",
  "Col5a3",
  "Col14a1"
)

lapply(as.list(list_genes_fibro_small), function(gene){
  umap_gene(sce_str_mmus, gene)
})

sce_fibro <- sce_str_mmus[rownames(sce_str_mmus) %in% c(
  list_genes_fibro_small
),]

# make a df for visualisation
logcounts_df_fibro <- logcounts(sce_fibro) %>%
  as.matrix() %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "genes") %>%
  tidyr::pivot_longer(
    cols = colnames(sce_fibro),
    names_to = "cell",
    values_to = "logcounts")

# add cluster info
logcounts_df_fibro$celltypes <- sce_fibro$celltypes[
  match(
    logcounts_df_fibro$cell,
    rownames(colData(sce_fibro)))]

logcounts_df_fibro$genes <- factor(
  logcounts_df_fibro$genes, 
  levels = c(list_genes_fibro_small))

vioplot_fibro <- logcounts_df_fibro %>%
  ggplot2::ggplot(
    aes(
      x = celltypes,
      y = logcounts,
      color = celltypes))+
    ggplot2::geom_violin(
      scale = "width",
      draw_quantiles = c(0.25, 0.5, 0.75))+
    #ggplot2::scale_y_continuous(
    #  limits = c(0, 6),
    #  breaks = c(0, 5)
    #)+
    theme_all+
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90, 
        hjust = 0.5,
        vjust = 0,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face),
      strip.text.y.right = element_text(
        angle = 0, 
        hjust = 0,
        vjust = 0.5,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face
      ),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none")+
    ggplot2::facet_grid(rows = vars(genes))+
    ggplot2::scale_color_manual("Cell types", values = col_cts_str)

vioplot_fibro


#------------------------------------------------------------------------------
# Osteo/balanced (Neo1)

list_genes_osteo <- c(
  "Neo1",

  # positive markers
  "Sdc3", # almost perfect in MSCs
  "Omd",
  "Ptprz1",
  "Bglap",
  "Spp1",

  # partial within population
  "Mmp13",
  "Timp1",
  "Car3",
  "Bmp3",

  # not exclusive 
  #"Cfh",
  #"Fgfr1",
  #"Col1a1",
  "Runx2"
  #"Col8a1"
)

# umaps

lapply(as.list(list_genes_osteo), function(gene){
  umap_gene(sce_str_mmus, gene)
})

# violin

sce_osteo <- sce_str_mmus[rownames(sce_str_mmus) %in% c(
  list_genes_osteo
),]

# make a df for visualisation
logcounts_df_osteo <- logcounts(sce_osteo) %>%
  as.matrix() %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "genes") %>%
  tidyr::pivot_longer(
    cols = colnames(sce_osteo),
    names_to = "cell",
    values_to = "logcounts")

# add cluster info
logcounts_df_osteo$celltypes <- sce_osteo$celltypes[
  match(
    logcounts_df_osteo$cell,
    rownames(colData(sce_osteo)))]

logcounts_df_osteo$genes <- factor(
  logcounts_df_osteo$genes, 
  levels = c(list_genes_osteo))

vioplot_osteo <- logcounts_df_osteo %>%
  ggplot2::ggplot(
    aes(
      x = celltypes,
      y = logcounts,
      color = celltypes))+
    ggplot2::geom_violin(
      scale = "width",
      draw_quantiles = c(0.25, 0.5, 0.75))+
    #ggplot2::scale_y_continuous(
    #  limits = c(0, 6),
    #  breaks = c(0, 5)
    #)+
    theme_all+
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90, 
        hjust = 0.5,
        vjust = 0,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face),
      strip.text.y.right = element_text(
        angle = 0, 
        hjust = 0,
        vjust = 0.5,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face
      ),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none")+
    ggplot2::facet_grid(rows = vars(genes))+
    ggplot2::scale_color_manual("Cell types", values = col_cts_str)

vioplot_osteo

# I would like to go with
# - Sdc3 observed in osteolineage cells
# - Omd major osteo marker
# - Ptprz1 major osteo marker
# - Bglap major osteo marker
# - Spp1 or Mmp13 major osteo markers


# violin small

list_genes_osteo_small <- c(
  "Neo1",
  "Sdc3", 
  "Omd", 
  "Ptprz1", 
  "Bglap", 
  "Spp1" 
)

lapply(as.list(list_genes_osteo_small), function(gene){
  umap_gene(sce_str_mmus, gene)
})

sce_osteo <- sce_str_mmus[rownames(sce_str_mmus) %in% c(
  list_genes_osteo_small
),]

# make a df for visualisation
logcounts_df_osteo <- logcounts(sce_osteo) %>%
  as.matrix() %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "genes") %>%
  tidyr::pivot_longer(
    cols = colnames(sce_osteo),
    names_to = "cell",
    values_to = "logcounts")

# add cluster info
logcounts_df_osteo$celltypes <- sce_osteo$celltypes[
  match(
    logcounts_df_osteo$cell,
    rownames(colData(sce_osteo)))]

logcounts_df_osteo$genes <- factor(
  logcounts_df_osteo$genes, 
  levels = c(list_genes_osteo_small))

vioplot_osteo <- logcounts_df_osteo %>%
  ggplot2::ggplot(
    aes(
      x = celltypes,
      y = logcounts,
      color = celltypes))+
    ggplot2::geom_violin(
      scale = "width",
      draw_quantiles = c(0.25, 0.5, 0.75))+
    #ggplot2::scale_y_continuous(
    #  limits = c(0, 6),
    #  breaks = c(0, 5)
    #)+
    theme_all+
    ggplot2::theme(
      axis.text.x = element_text(
        angle = 90, 
        hjust = 0.5,
        vjust = 0,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face),
      strip.text.y.right = element_text(
        angle = 0, 
        hjust = 0,
        vjust = 0.5,
        color = axis_text_color,
        size = axis_text_size,
        face = axis_text_face
      ),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none")+
    ggplot2::facet_grid(rows = vars(genes))+
    ggplot2::scale_color_manual("Cell types", values = col_cts_str)

vioplot_osteo

