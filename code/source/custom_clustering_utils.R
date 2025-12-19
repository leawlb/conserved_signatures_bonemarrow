
# annotations, colors, and genes for the custom species-specific clusterings

# since these are so many different cell types, I'm sourcing this out

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# HSC

#------------------------------------------------------------------------------

factors_custom_hsc <- c(
    "HSC",
    "HSC/Early Mk",
    "Early MPP cells 1",
    "Early MPP cells 2",
    "MPP cells 1",
    "MPP cells 2",
    "Activated MPP",
    "Activated MPP cells 1",
    "Activated MPP cells 2",
    "Late MPP",
    "Early myeloid",

    "Early granu/mono",
    "Granu/mono progenitor",
    "Neutro progenitor",
    "Later granu/neutro",

    "Early lymphoid",
    "Lymphoid",
    "Late lymphoid",
    "Late/Cycling lymphoid",

    "Baso/Mast progenitor",
    "Mk/E/BM (more Mk and BM)",
    "Mk/E/BM (more Ery)",
    "Early Mk progenitor",
    "Mk/Ery (more Mk)",
    "Mk/Ery",
    "Early ery/baso/mast 1",
    "Early ery/baso/mast 2",
    "Early ery",
    "Early erythroid progenitor",
    "Ery progenitor",
    "Erythroid 1",
    "Erythroid 2",
    "Erythroid 3",
    "Late erythroid progenitor",
    "Latest erythroid progenitor",
    "Cycling",
    "Cycling (more M-phase)",
    "Cycling (more S-phase)",
    "Cycling MPP",
    "Cycling, mostly myeloid",
    "Cycling granu/mono",
    "Cycling ery",
    # separate cells that would be removed
    "Antigen-presenting",
    "Mono/antigen-presenting",
    "Late neutro",
    "Late immune, mixed",
    "Late neutro/ery, mixed",
    "Low quality",
    "Low quality, early",
    "Low quality, neutro",
    "Low quality, ery"
)

factors_custom_hsc_unif <- c(
    "HSC",
    "HSC/Early Mk",
    "Early MPP cells",
    "MPP cells",
    "Activated MPP",
    "Activated MPP cells",
    "Late MPP",
    "Early myeloid",

    "Early granu/mono",
    "Granu/mono progenitor",
    "Neutro progenitor",
    "Later granu/neutro",

    "Early lymphoid",
    "Lymphoid",
    "Late lymphoid",
    "Late/Cycling lymphoid",

    "Baso/Mast progenitor",
    "Mk/E/BM (more Mk and BM)",
    "Mk/E/BM (more Ery)",
    "Early Mk progenitor",
    "Mk/Ery (more Mk)",
    "Mk/Ery",
    "Early ery/baso/mast",
    "Early ery",
    "Early erythroid progenitor",
    "Ery progenitor",
    "Erythroid",
    "Late erythroid progenitor",
    "Latest erythroid progenitor",
    "Cycling",
    "Cycling (more M-phase)",
    "Cycling (more S-phase)",
    "Cycling MPP",
    "Cycling, mostly myeloid",
    "Cycling granu/mono",
    "Cycling ery",
    # separate cells that would be removed
    "Antigen-presenting",
    "Mono/antigen-presenting",
    "Late neutro",
    "Late immune, mixed",
    "Late neutro/ery, mixed",
    "Low quality"
)

#------------------------------------------------------------------------------
# colors

col_custom_hsc <- c(
  "HSC" = "#4a0805ff", # mmus, mspr, mcar
  "HSC/Early Mk" = "#49006bff", # mcas
  "Early MPP cells 1" = "#6f120fff", # mmus, mcas
  "Early MPP cells 2" = "#b51611ff", # mmus, mcas
  "MPP cells 1" = "#fa221bff", # mspr
  "MPP cells 2" = "#ee8581ff", # mspr
  "Activated MPP" = "#c8605aff", # mmus, mcar 
  "Activated MPP cells 1" = "#f49792ff", # mcas
  "Activated MPP cells 2" = "#e95e57ff", # mcas
  "Late MPP" = "#e5530fff", # mcas
  "Early myeloid" = "#9d5c1fff", # mspr
  "Early granu/mono" = "#c87121ff", # mcas
  "Granu/mono progenitor" = "#e28c02ff", # mmus, mcar
  "Neutro progenitor" = "#f1ac3dff", # mspr, mcar
  "Later granu/neutro" = "#fc8c23ff", # mcas
  "Early lymphoid" = "#ffd000ff", # mmus, mspr
  "Lymphoid" = "#ffe28aff", # mcas
  "Late lymphoid" = "#faf24fff", # mcar
  "Late/Cycling lymphoid" = "#ffff75ff", # mspr
  "Baso/Mast progenitor" = "#ab00fbff", # mmus, mcas, mspr, mcar
  "Mk/E/BM (more Mk and BM)" = "#6e468bff", # mcar
  "Mk/E/BM (more Ery)" = "#9d81fbff", # mcar
  "Early Mk progenitor" = "#8d0153ff", # mspr
  "Mk/Ery (more Mk)" = "#8e3268ff", # mmus 
  "Mk/Ery" = "#d857a2ff", # mcas
  "Early ery/baso/mast 1" = "#993a8eff", # mspr
  "Early ery/baso/mast 2" = "#d857c9ff", # mspr
  "Early ery" = "#c519bdff", # mmus 
  "Early erythroid progenitor" = "#f4538eff", # mcar
  "Ery progenitor" = "#f989b2ff", # mmus, mcas, mspr
  "Erythroid 1" = "#b32e5fff", # mcar
  "Erythroid 2" = "#ff9bc0ff", # mcar
  "Erythroid 3" = "#f4bde6ff", # mcar
  "Late erythroid progenitor" = "#c26cacff", # mcar
  "Latest erythroid progenitor" = "#6e468bff", # mcar
  "Cycling" = "#b5b5b5ff", # mcas
  "Cycling (more M-phase)" = "#a68482ff", # mmus 
  "Cycling (more S-phase)" = "#cbb9b9ff", # mmus 
  "Cycling MPP" = "#c2ac89ff", # mcar
  "Cycling, mostly myeloid" = "#98877dff", # mspr
  "Cycling granu/mono" = "#dfdcd6ff", # mcar
  "Cycling ery" = "#c3a4b0ff", # mspr
  # separate cells that would be removed
  "Antigen-presenting" = "#4c64fbff", # mcar
  "Mono/antigen-presenting" = "#93a0f7ff", # mspr
  "Late neutro" = "#6cfff5ff", # mmus 
  "Late immune, mixed" = "#37ded3ff", # mmus, mcas, mspr
  "Late neutro/ery, mixed" = "#03b0a5ff", # mcas
  "Low quality" = "#4adf6dff", # mcar
  "Low quality, early" = "#208c39ff", # mspr
  "Low quality, neutro" = "#75811bff", # mspr
  "Low quality, ery" = "#99b835ff" # mspr
)

col_custom_hsc_unif <- c(
  "HSC" = "#4a0805ff", # mmus, mspr, mcar
  "HSC/Early Mk" = "#49006bff", # mcas
  "Early MPP cells" = "#b51611ff", # mmus, mcas
  "MPP cells" = "#ee8581ff", # mspr
  "Activated MPP" = "#c8605aff", # mmus, mcar 
  "Activated MPP cells" = "#e95e57ff", # mcas
  "Late MPP" = "#e5530fff", # mcas
  "Early myeloid" = "#9d5c1fff", # mspr
  "Early granu/mono" = "#c87121ff", # mcas
  "Granu/mono progenitor" = "#e28c02ff", # mmus, mcar
  "Neutro progenitor" = "#f1ac3dff", # mspr, mcar
  "Later granu/neutro" = "#fc8c23ff", # mcas
  "Early lymphoid" = "#ffd000ff", # mmus, mspr
  "Lymphoid" = "#ffe28aff", # mcas
  "Late lymphoid" = "#faf24fff", # mcar
  "Late/Cycling lymphoid" = "#ffff75ff", # mspr
  "Baso/Mast progenitor" = "#ab00fbff", # mmus, mcas, mspr, mcar
  "Mk/E/BM (more Mk and BM)" = "#6e468bff", # mcar
  "Mk/E/BM (more Ery)" = "#9d81fbff", # mcar
  "Early Mk progenitor" = "#8d0153ff", # mspr
  "Mk/Ery (more Mk)" = "#8e3268ff", # mmus 
  "Mk/Ery" = "#d857a2ff", # mcas
  "Early ery/baso/mast" = "#993a8eff", # mspr
  "Early ery" = "#c519bdff", # mmus 
  "Early erythroid progenitor" = "#f4538eff", # mcar
  "Ery progenitor" = "#f989b2ff", # mmus, mcas, mspr
  "Erythroid" = "#b32e5fff", # mcar
  "Late erythroid progenitor" = "#c26cacff", # mcar
  "Latest erythroid progenitor" = "#6e468bff", # mcar
  "Cycling" = "#b5b5b5ff", # mcas
  "Cycling (more M-phase)" = "#a68482ff", # mmus 
  "Cycling (more S-phase)" = "#cbb9b9ff", # mmus 
  "Cycling MPP" = "#c2ac89ff", # mcar
  "Cycling, mostly myeloid" = "#98877dff", # mspr
  "Cycling granu/mono" = "#dfdcd6ff", # mcar
  "Cycling ery" = "#c3a4b0ff", # mspr
  # separate cells that would be removed
  "Antigen-presenting" = "#4c64fbff", # mcar
  "Mono/antigen-presenting" = "#93a0f7ff", # mspr
  "Late neutro" = "#6cfff5ff", # mmus 
  "Late immune, mixed" = "#37ded3ff", # mmus, mcas, mspr
  "Late neutro/ery, mixed" = "#03b0a5ff", # mcas
  "Low quality" = "#4adf6dff", # mcar
  "Low quality, early" = "#208c39ff", # mspr
  "Low quality, neutro" = "#75811bff", # mspr
  "Low quality, ery" = "#99b835ff" # mspr
)

#------------------------------------------------------------------------------

# hsc marker genes list long
# comments relate to the species-specific clustering
  marker_genes_list_hsc <- list(
    "HSC" = c(
      "Ly6a", # HSC-specific, spotty
      "Slamf1",
      "Aldh1a1", 
      "Alcam",
      "Gng11",
      "Mllt3",
      "Socs2",

      "Cdkn1c", # HSC-specific, works for all and I known them well
      "Procr",
      "Mmrn1",
      "Mecom",
      "Mpl",
      "Pdzk1ip1", 

      "Ltb", # mostly HSC, works for all, and I know them
      "Hlf", 
      "Ldhb",

      "Txnip", # redox functions, likely also expressed in other clusters
      "Pbx1",

      "Angpt1", # quite broad already, but highest in HSC and well-known mostly
      "Hoxa9",
      "Lmo2",
      "Meis1",
      "Msi2"
    ),
    "MPP" = c(
      "Adgrl4", 
      "Cd34", 
      "Sox4",
      "Cd48", 
      "Sell",
      "Cdk4"
    ),
    "Myeloid_general" = c(
      "Spi1", 
      "Cebpa",
      "Ctsg",
      "Mpo",
      "Prtn3",
      "F13a1"),

    "Neutro" = c(
      "Elane", 
      "Gfi1",
      "Cebpe", 
      "S100a8",
      "S100a9",
      "Wfdc21",
      "Chil1"),

    "Mono" = c(  
      "Ccr2",
      "Irf8",
      "Ly86",
      "Ms4a6c"),

    "Antigen_presenting_contamination" = c(
      "Cd74",
      "H2-Aa",
      "H2-Eb1",
      "H2-Ab1",
      "Ctss"),

    "Lymphoid_early" = c(
      "Tcf4",
      "Flt3", 
      "Satb1",
      "Dntt"),

    "Lymphoid_late" = c(
      "Ighm",
      "Il7r",
      "Rag1",
      "Vpreb1",
      "Vpreb3",
      "Cd79b", 
      "Igkc"),

    "Lymph_contamination" = c(
      "Ccl5", # NK
      "Cd3e", # TC
      "Cd8a",
      "Tcf7"),

    "Mast/Baso" = c(
      "Mcpt8", 
      "Hdc", 
      "Lmo4", 
      "Gzmb",
      "Egr1",
      "Cpa3", 
      "Ms4a2",
      "Prss34",
      "Fcer1a",
      "Gata2"#  also Mk, Eo Baso and HSC
    ),

    "Mk" = c(
      "Itga2b", 
      "Cavin2",
      "Pf4",
      "Clec4d",
      "Plek",
      "Fli1"),

    "Mk/Ery" = c(
      "Zfpm1",
      "Gata1", 
      "Trib2"),

    "Ery" = c(
      "Frat2",
      "Ass1",
      "Aqp11",
      "Klf1",
      "Blvrb",
      "Tspo2",
      "Epor", 
      "Gfi1b",
      "Gypa",
      "Hba-a1",
      "Hbb-bs"),

    "CC entry/activation" = c(
      "Lig1",
      "Mcm7",
      "Hells",
      "Mcm6"),

    "Cell cycle S phase" = c(
      "Aurkb",
      "Cdk1"
      ),

    "Cell cycle M phase" = c(
      "Cdc20",
      "Ccnb2",
      "Knstrn"),

    "Cell cycle general" = c(
      "Mki67",
      "Lockd",
      "Cenpa" 
      )
    )

# for figure s3
# names related to merged clustering
marker_genes_list_hsc_figs3 <- list(
    "HSC_by_spec" = c(
      "Slamf1", # generally weak
      "Ly6a", # strong in BL6
      "Procr", 
      "Pdzk1ip1",
      "Tgm2", 
      "Socs2", # strong in cast
      "Aldh1a1", # strong in spret
      "Alcam",
      "Mllt3", # strong in caroli
      "Gng11",
      "Cited2"
    ),
    "HSC_shared" = c(
      "Mpl",
      "Mmrn1",
      "Mecom",
      "Cdkn1c", 
      "Txnip"
    ),  
    "Stem" = c( # this includes expression in MPPs, albeit lower
      "Ltb",
      "Lmo2",
      "Meis1",
      "Hlf", 
      "Hoxa9",
      "Angpt1", 
      "Msi2"
    ),
    "MPP" = c( 
      "Adgrl4", 
      "Cd34", 
      "Sox4",    
      "Flt3",
      "Cd48", 
      "Sell",
      "Cdk4"
    ),
    "CC_entry_activation" = c(
      "Lig1",
      "Mcm7",
      "Hells",
      "Mcm6"
    ),
    "Myeloid_general" = c( # as negative markers for these populations
      "Spi1", 
      "Cebpa")
  )

#------------------------------------------------------------------------------
# best_marker_genes, a shorter list

best_marker_genes_hsc <- c(
  "Ly6a", # stem
  "Procr",
  "Cdkn1c",
  "Mmrn1",
  "Mpl",
  "Mecom",
  "Hlf",
  "Meis1",
  "Msi2",
  "Adgrl4", # mpp compartment
  "Cd34",
  "Spi1", # early myeloid
  "Cebpa",
  "Mpo",
  "Elane",
  "Gfi1", # neutro
  "S100a8",
  "Wfdc21",
  "Irf8", # mono
  "Ly86",
  "Ms4a6c",
  "Flt3", # b early
  "Satb1",
  "Dntt",
  "Il7r",
  "Vpreb1", # b later
  "Cd79b",
  "Gzmb", # baso mast
  "Cpa3",
  "Ms4a2",
  "Cavin2", # mk
  "Pf4", 
  "Itga2b",
  "Zfpm1",
  "Klf1", # ery earliest
  "Tspo2",
  "Epor", # ery later
  "Gypa",
  "Hba-a1",
  "Lig1", # cc entry
  "Mcm6",
  "Aurkb", # s phase
  "Cdk1",
  "Cdc20", # mphase
  "Knstrn",
  "Mki67", # cycling general
  "Ccl5", # contamination, Nk
  "Cd8a",
  "Cd74", # contamination, antigen presenting
  "H2-Aa"
)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# STROMAL

#------------------------------------------------------------------------------

# str all custom labels correct order for factoring

factors_custom_str <- c(
  "Adipo CAR", 
  "Adipo/Osteo CAR", 
  "Osteo", 
      
  "Balanced mesenchymal", 

  "Early chondro",
  "Late chondro and Osteo",
  "Late chondro", 

  "Articular/synovial-like",

  "Fibro/Chondro (more Fibro)",
  "Fibro", 
      
  "Capillary EC", 
  "Capillary EC 1", 
  "Capillary EC 2", 
  "Arteriolar/capillary EC", 
  "Arteriolar EC",

  "Transitory/Venous EC", 
  "Mix arteriolar/sinus EC", 
  "Mix ven/sinus/lymph EC", 
  "Sinusoidal EC", 

  "Lymphatic EC", 

  "Peri/SMC", 

  "Erythroid", 
  "Ery 1", 
  "Ery 2", 
  "Ery 3", 
  "Ery 4", 
  "Ery 5", 
  "Ery 6", 
  "Ery 7", 
  "Mature ery",

  "Mono/Macro 1", 
  "Mono/Macro 2", 
  "Mono/Macro 3", 
  "Mono/Macro 4", 
  "Mono/Macro 5", 
  
  "Early Neutro 1", 
  "Early Neutro 2", 

  "Neutro 1", 
  "Neutro 2", 
  "Neutro 3",
  "Neutro 4",

  "Mature Neutro 1", 
  "Mature Neutro 2", 
  "Mature Neutro 3", 

  "Nk/TC", 
  "Eosinophil", 
  "Eo/Baso/Mast", 
  "Baso/Mast", 

  "BC lineage 1",
  "BC lineage 2", 
  "BC lineage 3", 
  "BC lineage 4", 
  "BC lineage 5", 
  "BC lineage 6", 
  "BC lineage 7", 

  "Low qu adipo CAR",
  "Low qu osteo-chondro",

  "Mature immune mix 1", 
  "Mature immune mix 2", 
  "Antigen-presenting",
  "Antigen-presenting 1", 
  "Antigen-presenting 2", 

  "Mix lymphatic EC/immune", 
  "Mix skeletal muscle/immune", 
  "Skeletal muscle",
  "Low qu skeletal muscle" 
)

factors_custom_str_unif <- c(
  "Adipo CAR", 
  "Adipo/Osteo CAR", 
  "Osteo", 
      
  "Balanced mesenchymal", 

  "Early chondro",
  "Late chondro and Osteo",
  "Late chondro", 

  "Articular/synovial-like",

  "Fibro/Chondro (more Fibro)",
  "Fibro", 
      
  "Capillary EC", 
  "Arteriolar/capillary EC", 
  "Arteriolar EC",

  "Transitory/Venous EC", 
  "Mix arteriolar/sinus EC", 
  "Mix ven/sinus/lymph EC", 
  "Sinusoidal EC", 

  "Lymphatic EC", 

  "Peri/SMC", 

  "Erythroid", 
  "Ery", 
  "Mature ery",

  "Mono/Macro", 
  "Early Neutro", 
  "Neutro", 
  "Mature Neutro", 

  "Nk/TC", 
  "Eosinophil", 
  "Eo/Baso/Mast", 
  "Baso/Mast", 

  "BC lineage",

  "Low qu adipo CAR",
  "Low qu osteo-chondro",

  "Mature immune mix", 
  "Antigen-presenting",

  "Mix lymphatic EC/immune", 
  "Mix skeletal muscle/immune", 
  "Skeletal muscle",
  "Low qu skeletal muscle" 
)

#------------------------------------------------------------------------------
# corresponding colors
col_custom_str <- c(
  "Adipo CAR" = "#2c8284ff", # mmus, mspr 
  "Adipo/Osteo CAR" = "#35abadff", # mcas, mcar
  "Osteo" = "#1bfcdeff", # mmus, mcas, mspr
    
  "Balanced mesenchymal" = "#64fc71ff", # mmus, mcas, mspr

  "Fibro" = "#557113ff",  # mmus, mcas, mcar
  "Fibro/Chondro (more Fibro)" = "#7bab0aff", # mspr
  "Early chondro" = "#c1fda3ff",  # mspr, mcar
  "Late chondro and Osteo" = "#80f000ff",  # mcar
  "Late chondro" = "#55a646ff", # mmus, mcas

  "Articular/synovial-like" = "#344b02ff",  # mspr, mcar

  "Capillary EC" = "#b8c1fbff", # mspr
  "Capillary EC 1" = "#b4e3fcff", # mmus
  "Capillary EC 2" = "#43bbfbff", # mmus
  "Arteriolar/capillary EC" = "#70abfeff", # mcas, mcar
  "Arteriolar EC" = "#2a66a2ff",  # mmus, mspr

  "Transitory/Venous EC" = "#6b44d6ff", # mmus, mspr
  "Mix arteriolar/sinus EC" = "#4c3f6eff", # mcar
  "Mix ven/sinus/lymph EC" = "#906eebff", # mcar
  "Sinusoidal EC" = "#352957ff", # mmus, mcas, mspr

  "Lymphatic EC" = "#4f7389ff", # mspr

  "Peri/SMC" = "#cc00ffff", # mmus, mcas, mspr, mcar

  "Erythroid" = "#ff005dff", # mmus, mcas
  "Ery 1" = "#e2407cff", # mcar
  "Ery 2" = "#ff92baff", # mcar
  "Ery 3" = "#a2003bff", # mcar
  "Ery 4" = "#610024ff", # mcar
  "Ery 5" = "#983056ff", # mcar
  "Ery 6" = "#ad5374ff", # mcar
  "Ery 7" = "#b87f94ff", # mcar
  "Mature ery" = "#b04747ff", # mcar

  "Mono/Macro 1" = "#ff5900ff", # mcas
  "Mono/Macro 2" = "#ffa97aff", # mcas
  "Mono/Macro 3" = "#bd754dff", # mcas
  "Mono/Macro 4" = "#c4581fff", # mcas
  "Mono/Macro 5" = "#8f3606ff", # mcas

  "Early Neutro 1" = "#ff8c00ff", # mcar
  "Early Neutro 2" = "#cf7f1dff", # mcar

  "Neutro 1" = "#ce8f42ff", # mcas, mcar
  "Neutro 2" = "#fdc98aff", # mcas, mcar
  "Neutro 3" = "#ceb493ff", # mcar
  "Neutro 4" = "#a6865fff", # mcar

  "Mature Neutro 1" = "#8e775aff", # mcar
  "Mature Neutro 2" = "#94662eff", # mcar
  "Mature Neutro 3" = "#78470bff", # mcar

  "Eosinophil" = "#ae2c42ff", # mcar
  "Eo/Baso/Mast" = "#aa2e32ff", # mspr
  "Baso/Mast"  = "#ff0008ff", # mcas, mcar

  "Nk/TC" = "#d0716cff", # mcas

  "BC lineage 1" = "#e3c036ff", # mcas
  "BC lineage 2" = "#ffe479ff", # mcas
  "BC lineage 3" = "#c3af62ff", # mcas
  "BC lineage 4" = "#958542ff", # mcas
  "BC lineage 5" = "#deec1cff", # mcas
  "BC lineage 6" = "#d9e07bff", # mcas
  "BC lineage 7" = "#9ea534ff", # mcas

  "Low qu adipo CAR" = "#acdfe0ff",  # mmus, mspr
  "Low qu osteo-chondro" = "#9ab8a3ff",  # mspr

  "Mature immune mix 1" = "#afa683ff", # mcas, mcar
  "Mature immune mix 2" = "#7b7660ff", # mcas, mcar
  "Antigen-presenting" = "#8b7467ff",  # mcar
  "Antigen-presenting 1" = "#655348ff", # mcas
  "Antigen-presenting 2" = "#4f4139ff", # mcas

  "Mix lymphatic EC/immune" = "#716757ff", # mmus
  "Mix skeletal muscle/immune" = "grey50", # mspr
  "Skeletal muscle" = "grey60",  # mcas
  "Low qu skeletal muscle" =  "grey80" # mmus
)

col_custom_str_unif <- c(
  "Adipo CAR" = "#2c8284ff", # mmus, mspr 
  "Adipo/Osteo CAR" = "#35abadff", # mcas, mcar
  "Osteo" = "#1bfcdeff", # mmus, mcas, mspr
    
  "Balanced mesenchymal" = "#64fc71ff", # mmus, mcas, mspr

  "Fibro" = "#557113ff",  # mmus, mcas, mcar
  "Fibro/Chondro (more Fibro)" = "#7bab0aff", # mspr
  "Early chondro" = "#c1fda3ff",  # mspr, mcar
  "Late chondro and Osteo" = "#80f000ff",  # mcar
  "Late chondro" = "#55a646ff", # mmus, mcas

  "Articular/synovial-like" = "#344b02ff",  # mspr, mcar

  "Capillary EC" = "#43bbfbff", # mspr
  "Arteriolar/capillary EC" = "#70abfeff", # mcas, mcar
  "Arteriolar EC" = "#2a66a2ff",  # mmus, mspr

  "Transitory/Venous EC" = "#6b44d6ff", # mmus, mspr
  "Mix arteriolar/sinus EC" = "#4c3f6eff", # mcar
  "Mix ven/sinus/lymph EC" = "#906eebff", # mcar
  "Sinusoidal EC" = "#352957ff", # mmus, mcas, mspr

  "Lymphatic EC" = "#4f7389ff", # mspr

  "Peri/SMC" = "#cc00ffff", # mmus, mcas, mspr, mcar

  "Erythroid" = "#ff005dff", # mmus, mcas
  "Ery" = "#ff92baff", # mcar
  "Mature ery" = "#b04747ff", # mcar

  "Mono/Macro" = "#ff5900ff", # mcas
  "Early Neutro" = "#cf7f1dff", # mcar
  "Neutro" = "#fdc98aff", # mcas, mcar
  "Mature neutro" = "#8e775aff", # mcar


  "Eosinophil" = "#ae2c42ff", # mcar
  "Eo/Baso/Mast" = "#aa2e32ff", # mspr
  "Baso/Mast"  = "#ff0008ff", # mcas, mcar

  "Nk/TC" = "#d0716cff", # mcas

  "BC lineage" = "#deec1cff", # mcas

  "Low qu adipo CAR" = "#acdfe0ff",  # mmus, mspr
  "Low qu osteo-chondro" = "#9ab8a3ff",  # mspr

  "Mature immune mix" = "#7b7660ff", # mcas, mcar
  "Antigen-presenting" = "#8b7467ff",  # mcar
  "Antigen-presenting" = "#4f4139ff", # mcas

  "Mix lymphatic EC/immune" = "#716757ff", # mmus
  "Mix skeletal muscle/immune" = "grey50", # mspr
  "Skeletal muscle" = "grey60",  # mcas
  "Low qu skeletal muscle" =  "grey80" # mmus
)

#------------------------------------------------------------------------------

# str genes in the correct order

marker_genes_list_str <- list(
    "CAR" = c( # classic niche 
      "Cxcl14", 
      "Ebf3", # niche cells https://doi.org/10.1101/gad.314013.118
      "Cxcl12", # also arteriolar EC
      "Kitl" # also arteriolar EC
      ),
    "Adipo" = c(
      "Adipoq",
      "Hp",
      "Lepr", 
      "Lpl",
      "Cebpa",
      "Pappa"  # adipocyte marker in MALP https://doi.org/10.7554/eLife.54695
      ),
    "Osteo" = c( # from early/CAR to more osteocyte-ish (approx)
      "Limch1", # Osteo-CAR https://doi.org/10.1016/j.jbc.2024.107158
      "Angpt4", # Osteo-CAR https://doi.org/10.1016/j.jbc.2024.107158
      "Wif1", # Osteo-CAR https://doi.org/10.1016/j.jbc.2024.107158
      "Spp1", # pre-osteoblast https://doi.org/10.1016/j.jbc.2024.107158
      "Mmp13", # pre-osteoblast https://doi.org/10.1016/j.jbc.2024.107158
      "Alpl", # Osteo-CAR https://doi.org/10.1038/s41556-019-0439-6
      "Slit2", # pre-osteoblast https://doi.org/10.1016/j.jbc.2024.107158
      "Sp7", # Osteo-CAR and osteoblast https://doi.org/10.1038/s41556-019-0439-6
      "Pth1r", # osteo ST1
      "Bglap", # osteo https://doi.org/10.1016/j.jbc.2024.107158
      "Bglap2", # osteoblasts https://doi.org/10.7554/eLife.54695
      "Bglap3", # 
      "Mepe", # osteocyte https://doi.org/10.1016/j.jbc.2024.107158
      "Ibsp", # osteoblasts https://doi.org/10.7554/eLife.54695
      "Phex" # osteocyte https://doi.org/10.1016/j.jbc.2024.107158
      ),
    "Chondro/Osteo/fibro/early mixed" = c(
      "Tnc", # osteo-X https://doi.org/10.1016/j.jbc.2024.107158
      "Ptgis", # osteo ST1
      "Postn", # osteo-X https://doi.org/10.1016/j.jbc.2024.107158
      "Ptn", # osteo-X https://doi.org/10.1016/j.jbc.2024.107158
      "Comp", # fibro-chondro ST1
      "Ogn", # in branch point cekk type to osteogenic state https://doi.org/10.1186/s13073-023-01224-0
      "Aspn", # osteo-X https://doi.org/10.1016/j.jbc.2024.107158
      "Gsn", # several undefined MSC populations 
      "Apod", # several undefined MSC populations 
      "Tnn" # osteo-X https://doi.org/10.1016/j.jbc.2024.107158
      ),
    "Chondro early/other/mixed" = c(
      "Fmod", # chondro ST1
      "Chad", # chondro ST1
      "Prelp", # chondro ST1
      "Mgp", # chondro ST1
      "Crispld1",  # art ch https://doi.org/10.1016/j.jbc.2024.107158
      "Cytl1" # art ch https://doi.org/10.1016/j.jbc.2024.107158
      ),
    "Articular cartilage/Synovial fibro" = c(  # these genes can be articular chondros or also early chondros or also early synovial fibros.
      "Htra4", # see mspr
      "Prg4", # art ch and syn fibro https://doi.org/10.1016/j.jbc.2024.107158
      "Creb5" # see mspr
      ),   
    "Chondro late" = c(  
      "Sox9", # chondro https://doi.org/10.1016/j.jbc.2024.107158
      "Snorc",  # chondro https://doi.org/10.1016/j.jbc.2024.107158
      "Acan"
      ),
    "Early Fibro mixed" = c(
      "Dcn", # fibro ST1
      "Lum", # fibro ST1
      "Dpt", # fibro https://doi.org/10.1016/j.jbc.2024.107158
      "Clec3b", # fibro https://doi.org/10.1016/j.jbc.2024.107158
      "Cd34", # also arteriolar EC
      "Thy1", # fibro, msc, related to osteogenesis https://doi.org/10.1096/fj.201701379R
      "Pdgfra", # fibro ST1
      "Pdgfrb" # 
      ),
    "Pan EC" = c(
      "Cdh5",
      "Pecam1"
      ),
    "EC other" = c( # COULD be interesting but mostly similarly expressed
      "Emcn",  # higher in H-type vessels
      "Kdr", # venous/sinusoid/unclear
      "Cd36"
      ),
   "Capillary" = c(
      "Car4",
      "Timp4",
      "Rgcc"
      ),
    "Arteriolar" = c(
      "Efnb2", 
      "Hey1", # have also seen as H type
      "Sox17" # have also seen as H type
      ),
    "Transitory/venous/fenestrated/other" = c(
      "Vwf", # venous or H type?
      "Ehd4", # angiogenesis or transitory?
      "Plvap" # fenestration
    ),
    "Sinusoidal/venous expressed" = c(
      "Selp",  # expressed in sinusoid
      "Icam1", # expressed in sinsuoid
      "Vcam1",# expressed in sinsuoid
      "Nr2f2", # venous identity
      "Ackr1", # expressed in sinusoid
      "Lrg1", # expressed in sinsuoid
      "Sele"  # expressed in sinusoid
    ),
    "Sinusoidal classic" = c(
      "Stab2", # classic sinusoid
      "Stab1", # classic sinusoid
      "Flt4", # classic sinusoid 
      "Dnase1l3", # sinusoid https://doi.org/10.1038/s41467-018-04726-3
      "Tfpi", # classic sinusoid
      "Sirpa", # sinusoid https://doi.org/10.1038/s41467-018-04726-3
      "C1qtnf1" # sinusoid https://doi.org/10.1038/s41467-018-04726-3
      ),
    "Lymphatic EC" = c(
      "Lyve1", # also sinusoid
      "Ccl21a"
      ),
    "Peri/SMC" = c(
      "Rgs5",
      "Acta2",
      "Tagln",
      "Steap4"
      ),
    "Contamination ery early/middle" = c(
      "Klf1",
      "Hmbs",
      "Hba-a1",
      "Hbb-bs",
      "Gypa"
      ),
    "Contamination ery late" = c(
      "Ermap",
      "Kel",
      "Sox6",
      "Slc4a1"
      ),
    "Contamination myeloid early" = c(
      "Mpo",
      "Elane",
      "Prtn3"
      ),
    "Contamination mono/macro" = c(
      "Irf8",
      "Ccr2",
      "F13a1",
      "Mpeg1"
      ),
    "Contamination neutro middle" = c(
      "S100a8",
      "S100a9",
      "Wfdc21",
      "Trem1",
      "Cebpe",
      "Camp",
      "Ngp",
      "Retnlg"
      ),
    "Contamination neutro late" = c(
      "Thbs1",
      "Ccl6"
      ),
    "Contamination eosinophil" = c(
      "Ear1",
      "Ear6",
      "Prg3"
    ),
    "Contamination baso mast" = c(
      "Prss34",
      "Gata2",
      "Hdc",
      "Cpa3"
      ),
    "Contamination NK lineage" = c(
      "Ccl5",
      "Nkg7"
      ), 
    "Contamination TC lineage" = c(
      "Cd8a",
      "Cd4",
      "Cd3e",
      "Cd3d",
      "Trbc1"
      ),
    "Contamination lymph early" = c(
      "Sox4",
      "Satb1",
      "Dntt",
      "Il7r"
      ),
    "Contamination BC lineage" = c(
      "Vpreb1",
      "Vpreb3",
      "Cd79a",
      "Cd79b",
      "Igha",
      "Jchain",
      "Igkc"
      ), 
    "Contamination Antigen-presenting" = c(
      "H2-Aa",
      "H2-Ab1",
      "Cd74",
      "Ctss"
      ),
    "Contamination DC" = c(
      "Siglech",
      "Cd7",
      "Runx2" # also osteo
      ), 
    "Contamination skeletal muscle" = c(
      "Acta1",
      "Tnnt3",
      "Tnnc2",
      "Tpm2"
    ),
    "CC entry/activation" = c(
      "Lig1",
      "Mcm7",
      "Hells",
      "Mcm6"
      ),
    "Cell cycle S phase" = c(
      "Aurkb",
      "Cdk1"
      ),
    "Cell cycle M phase" = c(
      "Cdc20",
      "Ccnb2",
      "Knstrn"
      ),
    "Cell cycle general" = c(
      "Mki67",
      "Lockd",
      "Cenpa" 
      )

  )

# sorted by literature but also by expression in our dataset
marker_genes_list_str_figs3 <- list(
    "Niche" = c( # classic niche 
      "Cxcl14", 
      "Ebf3", 
      "Cxcl12", 
      "Kitl" 
      ),
    "Adipo_CAR" = c(
      "Adipoq",
      "Lepr", 
      "Lpl",
      "Pappa" 
      ),
    "Osteo_CAR" = c( 
      "Limch1", 
      "Wif1", 
      "Alpl" 
    ),
    "Earlier_Osteo" = c(
      "Spp1", 
      "Mmp13",
      "Bglap", 
      "Ibsp" 
      ),
    "Early_mixed_expression" = c( 
      # These are genes that are expressed in some sort of multi-lineage
      # cells or progenitor/MSC-like populations.
      # Although they may be already strongly related to one lineage,
      # expression in early multi-lineage cells has been observed.
      # Mostly excluding adipo-lineage gene expression.
      # Ordered rom biased towards osteo to chondro to fibro-ish but excluding
      # the lineage exclusive genes AS OBSERVED IN OUT DATASET.
      "Tnc", # osteo-biased expression
      "Runx2",
      "Sp7", # osteo-chondro-biased epxression
      "Pth1r",
      "Ptgis", # osteo-fibro-biased expression
      "Postn",
      "Aspn", 
      "Lum",  # osteo-fibro-and-chondro-biased expression
      "Dcn",
      "Gsn", 
      "Apod", # fibro-chondro-biased
      "Comp", # chondro-biased expression
      "Prelp", 
      "Mgp"     
    ), 
    "Chondro" = c(
      "Fmod",
      "Chad", 
      "Crispld1",
      "Sox9",
      "Acan"
      ),
    "Articular_cartilage_synovial_fibro" = c(  
      "Htra4",
      "Prg4",
      "Creb5" 
      ),   
    "More_fibro" = c(
      "Dpt",
      "Clec3b", 
      "Cd34", 
      "Thy1"
      )
  )

# shorter list for str marker genes
best_marker_genes_str <- c(
  "Cxcl12", # niche
  "Kitl",
  "Adipoq", # adipo
  "Pappa",
  "Limch1", # osteo-CAR
  "Angpt4",
  "Mmp13", # osteo
  "Alpl",
  "Bglap",
  "Postn", # balanced/mixed
  "Aspn", 
  "Apod", # Mixed/Fibro
  "Gsn", 
  "Clec3b", # Fibro
  "Thy1",
  "Cd34",
  "Pdgfra",
  "Mgp", # Mixed/chondro
  "Fmod", 
  "Chad", # mid chondro
  "Sox9", # late chondro
  "Snorc", 
  "Prg4", # synovial/articular
  "Creb5", 
  "Cdh5", # pan-endo 
  "Pecam1", 
  "Car4", # capillaty
  "Timp4",
  "Rgcc",
  "Efnb2", # art
  "Hey1",
  "Sox17",
  "Vwf", # transitory/venous
  "Ehd4",
  "Plvap",
  "Selp", # transitory sinusoidal/venous
  "Nr2f2",
  "Ackr1", 
  "Stab2", # pure sinusoidal
  "Stab1",
  "Flt4",
  "Lyve1", # lymphatic
  "Ccl21a",
  "Rgs5", # peri SMC
  "Acta2",
  "Tagln",
  "Klf1", # ery
  "Hba-a1",
  "Ermap",
  "Kel",
  "Elane", # early-ish myeloid
  "Prtn3",
  "Ccr2", # mono
  "Mpeg1", # macro
  "S100a8", # early-ish neutro
  "s100a9",
  "Cebpe", # mid-ish neutro
  "Ngp",
  "Rentlg", # late neutro
  "Ccl6", # late neutro and eo
  "Ear1", # eo
  "Prg3",
  "Prss34",
  "Cpa3",
  "Nkg7", # Nk
  "Ccl5",
  "Cd8a", # TC
  "Cd4",
  "Dntt", # early lymph
  "Vpreb1", # mid-ish lymph
  "Vpreb3", 
  "Cd79a", 
  "Igha", # later lymph
  "Jchain",
  "H2-Aa", # antige,-presenting
  "H2-Ab1",
  "Cd74",
  "Siglech", # other immune
  "Acta1", # skeletal muscle
  "Tnnt3",
  "Tnnc2",
  "Lig1", # CC entry
  "Hells",
  "Aurkb", # S phase
  "Cdk1",
  "Ccnb2", # M phase
  "Knstrn",
  "Mki67",
  "Lockd",
  "Cenpa" 

)

# shortest list for str marker genes
minimal_marker_genes_str <- c(
  "Cxcl12", # niche
  "Kitl",
  "Adipoq", # adipo
  "Limch1", # osteo-CAR
  "Mmp13", # osteo
  "Bglap", 
  "Alpl",
  "Postn", # balanced/mixed
  "Mgp", # Mixed/chondro
  "Fmod", 
  "Chad", # mid chondro
  "Sox9", # late chondro
  "Snorc", 
  "Prg4", # synovial/articular
  "Creb5", 
  "Aspn", 
  "Gsn", # Mixed/Fibro
  "Clec3b", # Fibro
  "Cd34",
  "Pdgfra",
  "Cdh5", # pan-endo 
  "Pecam1", 
  "Car4", # capillaty
  "Timp4",
  "Rgcc",
  "Efnb2", # art
  "Hey1",
  "Sox17",
  "Vwf", # transitory/venous
  "Plvap",
  "Selp", # transitory sinusoidal/venous
  "Nr2f2",
  "Ackr1", 
  "Stab2", # pure sinusoidal
  "Flt4",
  "Lyve1", # lymphatic
  "Ccl21a",
  "Rgs5", # peri SMC
  "Acta2",
  "Klf1", # ery
  "Kel",
  "Elane", # early-ish myeloid
  "Prtn3",
  "Ccr2", # mono
  "Mpeg1", # macro
  "S100a8", # early-ish neutro
  "Cebpe", # mid-ish neutro
  "Rentlg", # late neutro
  "Ear1", # eo
  "Prss34",
  "Cpa3",
  "Nkg7", # Nk
  "Ccl5",
  "Cd8a", # TC
  "Cd4",
  "Dntt", # early lymph
  "Vpreb1", # mid-ish lymph
  "Cd79a", 
  "Igha", # later lymph
  "Jchain",
  "H2-Aa", # antige,-presenting
  "Cd74",
  "Siglech", # other immune
  "Acta1", # skeletal muscle
  "Tnnc2",
  "Lig1", # CC entry
  "Mki67", # general CC 
  "Lockd"
)
