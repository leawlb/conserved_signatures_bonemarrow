library(Seurat, quietly = TRUE)
library(ggplot2, quietly = TRUE)
set.seed(37)

source("../../source/sce_functions_reclustering.R")

data <- readRDS(snakemake@input[["data_input"]])
data.updated <- UpdateSeuratObject(object = data)  # available data is v3 Seurat
print(data)

# Optimize resolution using markers conserved across primate species
cons_markers <- readRDS(snakemake@input[["cons_markers"]])
seu_mark_cor <- unique(unlist(cons_markers)) # set of all conserved marker genes 

# to choose optimal cluster resolution later
resolution_vec <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)

# choosing resolution for reclustering with conserved markers
all_scores_markers <- list()

for(species in unique(data.updated$orig.ident)){
  
  seu <- subset(data.updated, subset = orig.ident == species)
  seu$cell_type <- as.factor(seu$subclass_label)
  
  seu_mark_cor_list <- lapply(X = as.list(resolution_vec),
                              FUN = standard_seu_pipeline,
                              seu = seu,
                              data_use = "raw_counts", # RNA assay
                              features = seu_mark_cor)
  names(seu_mark_cor_list) <- as.character(resolution_vec)
  
  scores <- lapply(X = seu_mark_cor_list,
                   FUN = calculate_scores)
  scores <- dplyr::bind_rows(scores, .id = "resolution")
  all_scores_markers[[species]] <- scores
}


saveRDS(all_scores_markers, file = snakemake@output[["all_scores_markers"]])



# choosing resolution for reclustering with signature genes

nDEGs <- readRDS(snakemake@input[["nDEGS"]])

# finding signature
conserved_signature <- lapply(names(cons_markers), function(x){
  markers <- cons_markers[[x]]
  ndeg <- nDEGs[[x]]
  return(markers[markers %in% ndeg])
})
names(conserved_signature) <- names(cons_markers)

seu_sign <- unique(unlist(conserved_signature)) # set of all conserved sign genes 


all_scores_signature <- list()

for(species in unique(data.updated$orig.ident)){
  
  seu <- subset(data.updated, subset = orig.ident == species)
  seu$cell_type <- as.factor(seu$subclass_label)
  
  seu_sign_list <- lapply(X = as.list(resolution_vec),
                              FUN = standard_seu_pipeline,
                              seu = seu,
                              data_use = "raw_counts", # RNA assay
                              features = seu_sign)
  names(seu_sign_list) <- as.character(resolution_vec)
  
  scores <- lapply(X = seu_sign_list,
                   FUN = calculate_scores)
  scores <- dplyr::bind_rows(scores, .id = "resolution")
  all_scores_signature[[species]] <- scores
}

saveRDS(all_scores_signature, file = snakemake@output[["all_scores_signature"]])
