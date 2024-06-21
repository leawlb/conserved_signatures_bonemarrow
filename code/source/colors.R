#-------------------------------------------------------------------------------

library(dittoSeq, quietly = TRUE)

#-------------------------------------------------------------------------------

# MARKERS
marker1 <- "dodgerblue3"
marker2 <- "darkgreen"

# used for ggplot
#scale_fill_continuous("% cells/cell type", 
#                      limits=c(0, 100), breaks=seq(0,100,by=20),
#                      low = "white", high = "blue")

# used for pheatmap (looks different!)
mycolors_to1 <- colorRampPalette(c("grey98", "blue"))(100)
names(mycolors_to1) <- c(1:100)/100

mycolors_to100 <- colorRampPalette(c("grey98", "blue"))(100)
names(mycolors_to100) <- c(1:100)


#-------------------------------------------------------------------------------

# SAMPLES/CLUSTER
col_num <- c(dittoColors()[2:35])
names(col_num) <- as.character(c(1:34))

col_alp <- c(dittoColors()[12:19])
names(col_alp) <- c("a", "b", "c", "d", "e", "f", "g")

#-------------------------------------------------------------------------------

# REFERENCE CELL TYPES
if(exists("colors_ref_path")){
  colors_ref_df <- read.csv(colors_ref_path, sep = ";", header = TRUE)
  col_cts_ref <- colors_ref_df$color
  names(col_cts_ref) <- colors_ref_df$celltype
}

# ALL OTHER COLORS 
if(exists("colors_path")){
  colors_df <- read.csv(colors_path, sep = ";", header = TRUE)
  
  col_age <- colors_df[colors_df$purpose == "Age_ID",]$color
  names(col_age) <- colors_df[colors_df$purpose == "Age_ID",]$level
  
  col_spc <- colors_df[colors_df$purpose == "Species_ID",]$color
  names(col_spc) <- colors_df[colors_df$purpose == "Species_ID",]$level
  
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
}
