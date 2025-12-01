
# annotations, colors, and genes for the custom species-specific 
# clusterings

# since these are so many different cell types, I'm sourcing this out

#------------------------------------------------------------------------------
# str all custom labels correct order for factoring

factors_custom_str <- c(
  "Adipo-CAR", 
  "Adipo/Osteo CAR", 
  "Osteo", 
      
  "Balanced mesenchymal", 

  "Early chondro",
  "Late chondro and Osteo",
  "Late chondro", 

  "Articular/synovial progenitor-like",

  "Fibro/Chondro (more Fibro)",
  "Fibro", 
      
  "Capillary EC", 
  "Capillary EC 1", 
  "Capillary EC 2", 
  "Arteriolar/capillary EC", 
  "Arteriolar EC",

  "Transitory/Venous EC", 
  "Mixed arteriolar/sinus EC", 
  "Mixed ven/sinus/lymph EC", 
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

  "Mono/Macro lineage 1", 
  "Mono/Macro lineage 2, cycling", 
  "Mono/Macro lineage 3", 
  "Mono/Macro lineage 4", 
  "Mono/Macro lineage 5", 
  
  "Early Neutro 1", 
  "Early Neutro 2", 

  "Neutro 1", 
  "Neutro 2", 
  "Neutro 3",
  "Neutro 4",

  "Mature neutro 1", 
  "Mature neutro 2", 
  "Mature neutro 3", 

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


  "Low quality adipo-CAR",
  "Low quality osteo-chondro",

  "Mature immune, mixed 1", 
  "Mature immune, mixed 2", 
  "Antigen-presenting",
  "Antigen-presenting 1", 
  "Antigen-presenting 2", 

  "Mixed lymphatic EC/mature immune", 
  "Immune/skeletal muscle, mixed", 
  "Skeletal muscle",
  "Low quality skeletal muscle" 
)

#------------------------------------------------------------------------------
# corresponding colors
col_custom_str <- c(
  "Adipo-CAR" = "#1c6c6eff", # mmus, mspr 
  "Adipo/Osteo CAR" = "#35abadff", # mcas, mcar
  "Osteo" = "#1bfcdeff", # mmus, mcas, mspr
    
  "Balanced mesenchymal" = "#64fc71ff", # mmus, mcas, mspr

  "Fibro" = "#557113ff",  # mmus, mcas, mcar
  "Fibro/Chondro (more Fibro)" = "#7bab0aff", # mspr
  "Early chondro" = "#c1fda3ff",  # mspr, mcar
  "Late chondro and Osteo" = "#80f000ff",  # mcar
  "Late chondro" = "#55a646ff", # mmus, mcas

  "Articular/synovial progenitor-like" = "#344b02ff",  # mspr, mcar

  "Capillary EC" = "#b8c1fbff", # mspr
  "Capillary EC 1" = "#b4e3fcff", # mmus
  "Capillary EC 2" = "#43bbfbff", # mmus
  "Arteriolar/capillary EC" = "#70abfeff", # mcas, mcar
  "Arteriolar EC" = "#2a66a2ff",  # mmus, mspr

  "Transitory/Venous EC" = "#6b44d6ff", # mmus, mspr
  "Mixed arteriolar/sinus EC" = "#4c3f6eff", # mcar
  "Mixed ven/sinus/lymph EC" = "#906eebff", # mcar
  "Sinusoidal EC" = "#22076eff", # mmus, mcas, mspr

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

  "Mono/Macro lineage 1" = "#ff5900ff", # mcas
  "Mono/Macro lineage 2, cycling" = "#ffa97aff", # mcas
  "Mono/Macro lineage 3" = "#bd754dff", # mcas
  "Mono/Macro lineage 4" = "#c4581fff", # mcas
  "Mono/Macro lineage 5" = "#8f3606ff", # mcas

  "Early Neutro 1" = "#ff8c00ff", # mcar
  "Early Neutro 2" = "#cf7f1dff", # mcar

  "Neutro 1" = "#ce8f42ff", # mcas, mcar
  "Neutro 2" = "#fdc98aff", # mcas, mcar
  "Neutro 3" = "#ceb493ff", # mcar
  "Neutro 4" = "#a6865fff", # mcar

  "Mature neutro 1" = "#8e775aff", # mcar
  "Mature neutro 2" = "#94662eff", # mcar
  "Mature neutro 3" = "#78470bff", # mcar

  "Nk/TC" = "#d0716cff", # mcas
  "Eosinophil" = "#ae2c42ff", # mcar
  "Eo/Baso/Mast" = "#aa2e32ff", # mspr
  "Baso/Mast"  = "#ff0008ff", # mcas, mcar

  "BC lineage 1" = "#e3c036ff", # mcas
  "BC lineage 2" = "#ffe479ff", # mcas
  "BC lineage 3" = "#c3af62ff", # mcas
  "BC lineage 4" = "#958542ff", # mcas
  "BC lineage 5" = "#deec1cff", # mcas
  "BC lineage 6" = "#d9e07bff", # mcas
  "BC lineage 7" = "#9ea534ff", # mcas

  "Low quality adipo-CAR" = "#acdfe0ff",  # mmus, mspr
  "Low quality osteo-chondro" = "#9ab8a3ff",  # mspr

  "Mature immune, mixed 1" = "#afa683ff", # mcas, mcar
  "Mature immune, mixed 2" = "#7b7660ff", # mcas, mcar
  "Antigen-presenting" = "#8b7467ff",  # mcar
  "Antigen-presenting 1" = "#655348ff", # mcas
  "Antigen-presenting 2" = "#4f4139ff", # mcas

  "Mixed lymphatic EC/mature immune" = "#716757ff", # mmus
  "Immune/skeletal muscle, mixed" = "#939393ff", # mspr
  "Skeletal muscle" = "#636363ff",  # mcas
  "Low quality skeletal muscle" =  "#acacacff" # mmus
)