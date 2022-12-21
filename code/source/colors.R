#-------------------------------------------------------------------------------

library(dittoSeq)

#-------------------------------------------------------------------------------

colors_celltypes <- read.csv(paste0(color_tables, "/colors_celltypes.csv"),
                             sep = ",")
colors_species <- read.csv(paste0(color_tables, "/colors_species.csv"),
                           sep = ",")
colors_annotations <- read.csv(paste0(color_tables, "/colors_annotations.csv"),
                               sep = ",")

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

# CELL TYPES
col_cts <- colors_celltypes[,3]
names(col_cts) <- colors_celltypes[,2]

# SPECIES
col_spc <- colors_species[,2]
names(col_spc) <- colors_species[,1]

# OTHER ANNOTATIONS
col_ann <- colors_annotations[,2]
names(col_ann) <- colors_annotations[,1]

#-------------------------------------------------------------------------------

col_list <- list(
  "Age" = col_ann,
  "Cluster" = col_num,
  "Species" = col_spc,
  #"Interaction Type" = col_itp,
  "Batch_exp" = col_alp,
  "Batch_seq" = col_num
)

#-------------------------------------------------------------------------------

# INTERACTION TYPES
#col_itp <- c(
#  "ECM" = "lemonchiffon2",
#  "ECM, Membrane" = "cornsilk3",
#  "ECM, Membrane, Secreted" = "honeydew3",
#  "ECM, Secreted" = "honeydew2",
#  "Membrane" = "darkslategray4",
#  "Membrane, Secreted" = "cadetblue3",
#  "Secreted" = "cadetblue2",
#  "Other" = "orange"
#)
