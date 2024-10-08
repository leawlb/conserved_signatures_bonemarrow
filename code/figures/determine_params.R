
# determine global params 

#### TEXT #### 

legend_title_face <- "plain"
legend_title_size <- 16
legend_title_color <- "black"

legend_text_face <- "plain"
legend_text_size <- 14
legend_text_color <- "black"

axis_title_face <- "plain"
axis_title_size <- 16
axis_title_color <- "black"

axis_text_face <- "plain"
axis_text_size <- 14
axis_text_color <- "black"
axis_text_size_small <- 10

plot_title_face <- "plain"
plot_title_size <- 16
plot_title_color <- "black"

#### OTHER THEME STUFF

strip_background_color <- "white"
strip_background_fill <- "white"
axis_ticks_color <- "black"


#### POINTS ####

umap_point_size <- 0.3
umap_point_alpha <- 1
umap_legend_point_size <- 3
umap_legend_point_alpha <- 1


max_size_dotplots <- 4



# specify one theme for all, customize elements per plot
theme_all <- ggplot2::theme_classic()+
  ggplot2::theme(
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

