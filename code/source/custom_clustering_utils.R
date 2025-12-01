
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
