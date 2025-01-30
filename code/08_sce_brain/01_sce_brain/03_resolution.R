library(Seurat, quietly = TRUE)
library(ggplot2, quietly = TRUE)

source("../../source/sce_functions_reclustering.R")

#data <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS")
data <- readRDS(snakemake@input[["data_input"]])
data.updated <- UpdateSeuratObject(object = data)  # available data is v3 Seurat
print(data)

# Optimize resolution using markers conserved across primate species
#core_markers <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/08_sce_brain/02_marker_conserved_primates.rds")
core_markers <- readRDS(snakemake@input[["core_markers"]])
seu_mark_cor <- unique(unlist(core_markers)) # set of marker genes 

# to choose optimal cluster resolution later
resolution_vec <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
                    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)

all_scores <- list()

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
  all_scores[[species]] <- scores
}

#saveRDS(all_scores,
#        "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/08_sce_brain/03_resolution_scores.rds")
saveRDS(all_scores, file = snakemake@output[["all_scores"]])


# ggplot(all_scores[["human"]],
#        aes(x = resolution,
#            y = value,
#            group = type)) +
#   facet_wrap(~type,
#                       scale = "free_y",
#                       ncol = 1) +
#   geom_line() +
#   theme_classic()
## human is 0.5

# ggplot(all_scores[["macaque"]],
#                 aes(x = resolution,
#                              y = value,
#                              group = type)) +
#   facet_wrap(~type,
#                       scale = "free_y",
#                       ncol = 1) +
#   geom_line() +
#   theme_classic()
# # macaque is 0.7

# ggplot(all_scores[["marmoset"]],
#                 aes(x = resolution,
#                              y = value,
#                              group = type)) +
#   facet_wrap(~type,
#                       scale = "free_y",
#                       ncol = 1) +
#   geom_line() +
#   theme_classic()
# # marmoset is 0.6

# ggplot(all_scores[["mouse"]],
#                 aes(x = resolution,
#                              y = value,
#                              group = type)) +
#   facet_wrap(~type,
#                       scale = "free_y",
#                       ncol = 1) +
#   geom_line() +
#   theme_classic()
# # mouse is 0.7
