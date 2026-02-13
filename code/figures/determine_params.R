
# determine global params 

#### TEXT #### 

legend_title_face <- "plain"
legend_title_size <- 20
legend_title_color <- "black"

legend_text_face <- "plain"
legend_text_size <- 20
legend_text_color <- "black"

axis_title_face <- "plain"
axis_title_size <- 20
axis_title_color <- "black"

axis_text_face <- "plain"
axis_text_size <- 20
axis_text_color <- "black"
axis_text_size_small <- 15

plot_title_face <- "plain"
plot_title_size <- 24
plot_title_color <- "black"

#### COLOR

strip_background_color <- "white"
strip_background_fill <- "white"
axis_ticks_color <- "black"

#### POINTS ####

umap_point_size <- 0.1
umap_point_alpha <- 1
umap_legend_point_size <- 3
umap_legend_point_alpha <- 1

max_size_dotplots <- 4

#### THEME ####

# specify one theme for all, customize elements per plot
theme_all <- ggplot2::theme_classic()+
  ggplot2::theme(
    plot.background = element_blank(),
    legend.background = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(
      size = axis_text_size,
      face = axis_text_face,
      color = axis_text_color),
    axis.title = element_text(
      size = axis_title_size,
      face = axis_title_face,
      color = axis_title_color),
    plot.title = element_text(
      size = plot_title_size,
      face = plot_title_face,
      color = plot_title_color),
    legend.text = element_text(
      size = legend_text_size,
      face = legend_text_face,
      color = legend_text_color),
    legend.title = element_text(
      size = legend_title_size,
      face = legend_title_face,
      color = legend_title_color),
    axis.ticks = element_line(
      color = axis_ticks_color))

theme_all_supp <- ggplot2::theme_classic()+
  ggplot2::theme(
    plot.background = element_blank(),
    legend.background = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(
      size = 16,
      face = axis_text_face,
      color = axis_text_color),
    axis.title = element_text(
      size = 16,
      face = axis_title_face,
      color = axis_title_color),
    plot.title = element_text(
      size = 16,
      face = plot_title_face,
      color = plot_title_color),
    legend.text = element_text(
      size = 16,
      face = legend_text_face,
      color = legend_text_color),
    legend.title = element_text(
      size = 16,
      face = legend_title_face,
      color = legend_title_color),
    axis.ticks = element_line(
      color = axis_ticks_color))



# determine some colors here ONLY for supplementary figures! so it
# fits best to this figure-specific source file.
# regarding the diffently reclustered other datasets for figure S5.
# it doesn't fit into source/custom_clustering_utils.R because that 
# is only for our own datas
# but the same or similar colors will be used.

col_anno_ts_bone_marrow <- c(
  "HSC" = "#4a0805ff",
  "HSC1 10xV1" = "#6d4545ff",
  "HSC2 10xV2" = "#d78585ff",
  "HSC3" = "#b41c14ff",
  "HSC, female-biased" = "#df7c77ff",
  "HSC, male-biased" = "#9a1f07ff",
  "Cycling Early 10xV1" =  "#958584ff",
  "Cycling 10xV2" =  "#dfdcd6ff",
  "Cycling" =  "#b5b5b5ff",
  "Early Mk prog" = "#4b003dff", 
  "Mk prog" = "#8e3268ff", 
  "Mk/Ery prog" = "#d857a2ff", 
  "Early Ery prog" = "#d195abff", 
  "Ery prog" = "#f4bde6ff",
  "Ery prog, female-biased" = "#f0c2f3ff",
  "Ery prog, male-biased" = "#946494ff",
  "Late Ery prog" = "#6e468bff",
  "Late ery prog, male-biased" = "#43047eff",
  "Baso/Eo/Mast/Granu prog" = "#ff0008ff", 
  "Early Myeloid" = "#774d23ff", 
  "Early Myeloid prog" = "#774d23ff", 
  "Myeloid prog" = "#ff9e02ff", 
  "Mono prog" = "#ff5900ff", 
  "Granu/Mono prog" = "#df8a02ff", 
  "Early granu prog, male-biased" = "#bb904aff", 
  "Early Granu/Mono prog" = "#937039ff", 
  "Mono/Dendritic prog" = "#ff5900ff",
  "Granu prog" = "#feb168ff",
  "Neutro prog, male-biased" = "#855f00ff",
  "Neutro prog" = "#fdda83ff",
  "More Granu prog" = "#f3cba5ff",
  "More Mono prog" = "#ca3504ff",
  "Lymphoid prog" = "#e1d337ff", 
  "Lymphoid MPP" = "#978c1aff",
  "Cycling Myeloid" = "#ba9687ff",
  "Cycling Myeloid prog" = "#ba9687ff",
  "batch, Smart-Seq2" = "#4f4cfdff",
  "Mk/Ery1 10xV1" = "#892e64ff",
  "Mk/Ery2 10xV2" = "#ff94d3ff",
  "Mk/Ery" = "#d857a2ff",
  "Antigen-presenting and lymphoid" = "#c3dd1eff"
)

col_anno_mus_tik_stromal <- c(
  "P1 adipogenesis1" = "#0f6f65ff",
  "P2 adipogenesis2" = "#1dcab9ff",
  "P3 osteo early" = "#95e9e9ff",
  "P4 osteo later" = "#74b2a8ff",
  "O1 osteo, Col16a" = "#56736fff",
  "O2 osteo-chondro" = "#92ff16ff",
  "O3 osteoblasts" = "#477e62ff",
  "V1 arterial" = "#4796e4ff",
  "V2 sinusoidal" = "#472a99ff",

  "Adipo/CAR1, pbs" = "#004133ff",
  "Adipo/CAR2, pbs"  = "#009f7dff",
  "Adipo/CAR3, non-pbs"  = "#00ffc3ff",
  "Osteo/CAR" = "#69cebdff",
  "Osteo-Chondro" = "#92ff16ff",
  "Fibro-Chondro, pbs" = "#a7de8cff",
  "Chondro" = "#55a646ff",
  "Osteo1" = "#1bf4fcff",
  "Osteo2" = "#76afb0ff",
  "Osteo3" = "#628182ff",
  "Art/Cap-EC" = "#4793ffff",
  "Ven EC1, non-pbs" = "#4a11e7ff",
  "Ven EC2" = "#8f6bf5ff",
  "Ven EC3, non-pbs" = "#b8add5ff"
)
