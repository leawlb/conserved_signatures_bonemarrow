
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
axis_text_size_small <- 14

plot_title_face <- "plain"
plot_title_size <- 24
plot_title_color <- "black"

#### OTHER THEME STUFF

strip_background_color <- "white"
strip_background_fill <- "white"
axis_ticks_color <- "black"


#### POINTS ####

umap_point_size <- 0.1
umap_point_alpha <- 1
umap_legend_point_size <- 3
umap_legend_point_alpha <- 1


max_size_dotplots <- 4



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

# TODO: put in colors later
colors_z_score <- c("blue3", "white", "red3")

col_cons <- c("conserved_signature" = "#FF4900",
              "conserved_markers" = "#FF9E1C",
              "mmusall_markers" = "#2DCA67",
              "ndges" = "#ea91adff",
              "random_features" = "grey80")

col_cons_long <- c(
  "conserved identity signature" = "#FF4900",
  "conserved markers" = "#FF9E1C",
  "BL6 markers" = "#2DCA67",
  "ndges" = "#ea91adff",
  "random genes" = "grey80",
  "species-specific marker genes" = "grey50")

col_spc_pub <- c(
  "BL6" = "black",
  "CAST" = "grey70",
  "SPRET" = "#B35100",
  "CAROLI" = "#FCC961"
)
