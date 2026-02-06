
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
  "HSC1 10xV1" = "#280604ff",
  "HSC2 10xV2" = "#a36363ff",
  "HSC3" = "#80201bff",
  "Cycling Early 10xV1" =  "#958584ff",
  "Cycling 10xV2" =  "#dfdcd6ff",
  "Cycling" =  "#b5b5b5ff",
  "Mk prog" = "#8e3268ff", 
  "Mk/Ery prog" = "#d857a2ff", 
  "Early Ery prog" = "#d195abff", 
  "Ery prog" = "#f4bde6ff",
  "Late Ery prog" = "#6e468bff",
  "Baso/Eo/Mast/Granu prog" = "#ff0008ff", 
  "Early Myeloid" = "#633c16ff", 
  "Early Myeloid prog" = "#633c16ff", 
  "Myeloid prog" = "#ff9e02ff", 
  "Mono prog" = "#ff5900ff", 
  "Granu/Mono prog" = "#df8a02ff", 
  "Early Granu/Mono prog" = "#937039ff", 
  "Mono/Dendritic prog" = "#ff5900ff",
  "Granu prog" = "#feb168ff",
  "More Granu prog" = "#f3cba5ff",
  "More Mono prog" = "#ca3504ff",
  "Lymphoid prog" = "#e1d337ff", 
  "Lymphoid MPP" = "#978c1aff",
  "Cycling Myeloid" = "#ba9687ff",
  "Cycling Myeloid prog" = "#ba9687ff",
  "batch, Smart-Seq2" = "#4f4cfdff",
  "Mk/Ery1 10xV1" = "#892e64ff",
  "Mk/Ery2 10xV2" = "#ff94d3ff",
  "Mk/Ery" = "#d857a2ff"
)
