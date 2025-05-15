# defining some color schemes here
# some of it based on metadata csv especially for specific cell types
# or similar

#-------------------------------------------------------------------------------

library(dittoSeq, quietly = TRUE)

#-------------------------------------------------------------------------------

# MARKERS
marker1 <- "dodgerblue3"
marker2 <- "darkgreen"

# used for pheatmap (looks different!)
mycolors_to1 <- grDevices::colorRampPalette(c("grey98", "blue"))(100)
names(mycolors_to1) <- c(1:100)/100

mycolors_to1_grey <- grDevices::colorRampPalette(c("grey85", "blue"))(100)
names(mycolors_to1_grey) <- c(1:100)/100

mycolors_to100 <- grDevices::colorRampPalette(c("grey98", "blue"))(100)
names(mycolors_to100) <- c(1:100)

mycolors_to_max <- grDevices::colorRampPalette(c("grey80", "#073A91"))(100)
names(mycolors_to_max) <- c(1:100)/100

# used for visualising QC (sum and detected)
color_vector_sum <- c("black", "darkorange4", "darkorange3", "darkorange1",
                      "orange", "orange", "gold1", "gold",
                      "lightgoldenrod", "lightgoldenrod",  "lightgoldenrod",
                      "lightgoldenrod", "lightgoldenrod", "lightgoldenrod")
color_vector_det <- c("black", "darkorange3", "orange", "lightgoldenrod")

col_mode <-c("marker" = "#EDE2D4",
             "conserved_marker" = "#EB9523",
             "conserved_signature" = "#DD5C00")
col_mode_log <-c("FALSE" = "#EDE2D4",
                 "TRUE" = "#DD5C00")

col_cons <- c("conserved_signature" = "#C8451D",
              "conserved_markers" = "#E48204",
              "mmusall_markers" = "#94B88E",
              "ndges" = "#4684FF",
              "random_features" = "#C3C3C3")

col_cons_long <- c(
  "conserved identity signature" = "#C8451D",
  "conserved markers" = "#E48204",
  "BL6 markers" = "#94B88E",
  "ndges" = "#4684FF",
  "random genes" = "#C3C3C3",
  "species-specific marker genes" = "grey50")

colors_z_score <- c("steelblue3", "white", "red3")

#-------------------------------------------------------------------------------

# SAMPLES/CLUSTER
col_num <- c(dittoSeq::dittoColors()[2:35])
names(col_num) <- as.character(c(1:34))

col_alp <- c(dittoSeq::dittoColors()[12:19])
names(col_alp) <- c("a", "b", "c", "d", "e", "f", "g")

#-------------------------------------------------------------------------------

# REFERENCE CELL TYPES
if(exists("colors_ref_path")){
  colors_ref_df <- utils::read.csv(colors_ref_path, sep = ";", header = TRUE)
  col_cts_ref <- colors_ref_df$color
  names(col_cts_ref) <- colors_ref_df$celltype
}

# ALL OTHER COLORS 
if(exists("colors_path")){
  colors_df <- utils::read.csv(colors_path, sep = ";", header = TRUE)
  
  col_age <- colors_df[colors_df$purpose == "Age_ID",]$color
  names(col_age) <- colors_df[colors_df$purpose == "Age_ID",]$level
  
  col_spc <- colors_df[colors_df$purpose == "Species_ID",]$color
  names(col_spc) <- colors_df[colors_df$purpose == "Species_ID",]$level
  
  col_spc_pub <- colors_df[colors_df$purpose == "Species_ID",]$color
  names(col_spc_pub) <- c("BL6", "CAST", "SPRET", "CAROLI")

  col_frc <- colors_df[colors_df$purpose == "Fraction_ID",]$color
  names(col_frc) <- colors_df[colors_df$purpose == "Fraction_ID",]$level
  
  col_asn <- colors_df[colors_df$purpose == "Assignment",]$color
  names(col_asn) <- colors_df[colors_df$purpose == "Assignment",]$level
  
  # cell types per fraction
  colors_df_temp <- colors_df[colors_df$purpose == "celltypes",]
  
  col_cts_hsc <- colors_df_temp[colors_df_temp$fraction == "hsc",]$color
  names(col_cts_hsc) <- colors_df_temp[colors_df_temp$fraction == "hsc",]$level
  
  col_cts_str <- colors_df_temp[colors_df_temp$fraction == "str",]$color
  names(col_cts_str) <- colors_df_temp[colors_df_temp$fraction == "str",]$level
  
  # categories per fraction
  colors_df_temp <- colors_df[colors_df$purpose == "category",]
  
  col_cat_hsc <- colors_df_temp[colors_df_temp$fraction == "hsc",]$color
  names(col_cat_hsc) <- colors_df_temp[colors_df_temp$fraction == "hsc",]$level
  
  col_cat_str <- colors_df_temp[colors_df_temp$fraction == "str",]$color
  names(col_cat_str) <- colors_df_temp[colors_df_temp$fraction == "str",]$level
  
  # antibodies per fraction
  colors_df_temp <- colors_df[colors_df$purpose == "Antibody_combination",]
  col_ab_hsc <- colors_df_temp[colors_df_temp$fraction == "hsc",]$color
  names(col_ab_hsc) <- colors_df_temp[colors_df_temp$fraction == "hsc",]$level
  
  col_ab_str <- colors_df_temp[colors_df_temp$fraction == "str",]$color
  names(col_ab_str) <- colors_df_temp[colors_df_temp$fraction == "str",]$level
  
}
