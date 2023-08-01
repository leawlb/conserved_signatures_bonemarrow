#-------------------------------------------------------------------------------

library(dittoSeq, quietly = TRUE)

#-------------------------------------------------------------------------------

# MARKERS
marker1 <- "dodgerblue3"
marker2 <- "darkgreen"

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
