#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

library(dittoSeq)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# general color scheme

# markers
marker1 <- "dodgerblue3"
marker2 <- "darkgreen"

# AGE
col_age <- c(
  "young" = "black",
  "old" = "grey62"
)

# SPECIES
col_spc <- c(
  "mmus" = "grey20",
  "mcas" = "lightsalmon4",
  "mcar" = "bisque4",
  "mspr" = "tan4"
)

# ANNOTATION
col_ann <- c(
  "emitter" = "darkslategray3",
  "receiver" = "indianred2"
)
col_lrs <- c(
  "ligand" = "darkslategray3",
  "receptor" = "indianred2"
)

# SAMPLES/CLUSTER
col_num <- c(dittoColors()[2:21])
names(col_num) <- as.character(c(1:20))

col_alp <- c(dittoColors()[12:19])
names(col_alp) <- c("a", "b", "c", "d", "e", "f", "g")

col_cls_man <- c(
  "1" = "grey30",
  "2" = "lightgoldenrod2",
  "3" = "#99CCFF",
  "4" = "grey80",
  "5" = "#FF6666",
  "6" = "#990000"
)

col_cls_man_anno <- c(
  "mixed young" = "grey30",
  "mixed old" = "grey80",
  "mixed ECs" = "#99CCFF",
  "old LT-HSCs" = "#990000",
  "old Mono/Dendr. progs" = "#FF6666",
  "mixed LMPPs" = "lightgoldenrod2"
)


# INTERACTION TYPES
col_itp <- c(
  "ECM" = "lemonchiffon2",
  "ECM, Membrane" = "cornsilk3",
  "ECM, Membrane, Secreted" = "honeydew3",
  "ECM, Secreted" = "honeydew2",
  "Membrane" = "darkslategray4",
  "Membrane, Secreted" = "cadetblue3",
  "Secreted" = "cadetblue2",
  "Other" = "orange"
)

#-------------------------------------------------------------------------------

# CELL TYPES

col_cts_all <- c(
  
  "Ng2-MSCs" = "darkslategray",
  "Fibro/Chondro progs" = "aquamarine4",
  "Arteriolar fibros" = "seagreen3",
  "Endosteal fibros" = "palegreen2",
  "Stromal fibros" = "darkolivegreen1",
  "Myofibroblasts" = "palegreen1",
  
  "Adipo-CARs" = "plum2",
  "Osteo-CARs" = "darkorchid1",
  
  "Arteriolar ECs" = "deepskyblue2",
  "Sinusoidal ECs" = "blue3",
  "Smooth muscle" = "lightblue1",
  
  "Chondrocytes" = "aquamarine",
  "Osteoblasts" = "turquoise1",
  
  "Schwann cells" = "lightcyan",
  
  "LT-HSCs" = "grey10",
  "ST-HSCs" = "grey30",
  "Multipotent prog, cycling" = "grey50",
  "ST-HSCs/MPPs" = "grey40",
  
  "LMPPs" = "wheat4",
  
  "Small pre-BCs" = "navajowhite3",
  "Large pre-BCs" = "navajowhite2",
  "Pro-BCs" = "antiquewhite1",
  "Immature B cells" = "lightyellow2",
  "Mature B cells" = "lightyellow1",
  "T cells/NK cells" = "snow3",
  
  "Common myeloid progs" = "sienna4",
  "Gran/Mono progs" = "sienna3",
  "Monocyte-primed progs" = "red3",  
  "Monocyte progs" = "red3",
  "Dendritic-primed progs" = "coral3",  
  "Neutrophil progs" = "darkorange2",
  "Baso/Eo/Mast progs" = "gold",
  
  "Monocytes" = "firebrick1",
  "Dendritic cells" = "lightsalmon",
  "Neutrophils" = "orange",
  
  "Ery/Mk progs" = "deeppink4",
  "Megakaryocyte progs" = "violetred1",
  "Erythrocyte progs" = "lightpink3",
  "Erythroid" = "pink"
)

col_cts_red <- c(

  "Ng2-MSCs" = "darkslategray",
  "Fibro/Chondro progs" = "turquoise4",
  "Arteriolar fibros" = "seagreen3",
  "Endosteal fibros" = "palegreen2",
  "Stromal fibros" = "darkolivegreen1",
  
  "Adipo-CARs" = "darkorchid1",
  "Osteo-CARs" = "plum1",
  
  "Arteriolar ECs" = "blue3",
  "Sinusoidal ECs" = "deepskyblue2",
  "Smooth muscle" = "lightblue1",
  
  "Osteoblasts" = "turquoise1",
  
  "LT-HSCs" = "grey10",
  "ST-HSCs/MPPs" = "grey40",
  
  "LMPPs" = "wheat4",
  
  "Common myeloid progs" = "sienna4",
  "Monocyte-primed progs" = "orangered3",  
  "Dendritic-primed progs" = "coral",  
  "Neutrophil progs" = "orange",
  "Baso/Eo/Mast progs" = "gold",
  
  "Ery/Mk progs" = "deeppink4",
  "Megakaryocyte progs" = "violetred1",
  "Erythrocyte progs" = "pink1"
)

col_rec_red <- c(
  "LT-HSCs" = "grey10",
  "ST-HSCs/MPPs" = "grey40",

  "LMPPs" = "wheat4",
  
  "Common myeloid progs" = "sienna4",
  "Monocyte-primed progs" = "orangered3",  
  "Dendritic-primed progs" = "coral",  
  "Neutrophil progs" = "orange",
  "Baso/Eo/Mast progs" = "gold",
  
  "Ery/Mk progs" = "deeppink4",
  "Megakaryocyte progs" = "violetred1",
  "Erythrocyte progs" = "pink1"
)

col_emi_red <- c(
  "Ng2-MSCs" = "darkslategray",
  "Fibro/Chondro progs" = "turquoise4",
  "Arteriolar fibros" = "seagreen3",
  "Endosteal fibros" = "palegreen2",
  "Stromal fibros" = "darkolivegreen1",
  
  "Adipo-CARs" = "darkorchid1",
  "Osteo-CARs" = "plum1",
  
  "Arteriolar ECs" = "blue3",
  "Sinusoidal ECs" = "deepskyblue2",
  "Smooth muscle" = "lightblue1",
  
  "Osteoblasts" = "turquoise1"
)

#-------------------------------------------------------------------------------

col_list <- list(
  "Age" = col_age,
  "Cluster" = col_num,
  "Species" = col_spc,
  "Interaction Type" = col_itp,
  "Celltypes" = col_cts_red,
  "Emitters" = col_rec_red,
  "Receivers" = col_emi_red,
  "Batch_exp" = col_alp,
  "Batch_seq" = col_num
)
