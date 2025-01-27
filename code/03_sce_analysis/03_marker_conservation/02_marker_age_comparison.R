## separate old and young, look for markers
## compare markers from old only, young only, and pooled data

library(Seurat)
library(tidyverse)

###### Hematapoetic
#data_hsc <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10")
#data_hsc <- readRDS("/omics/odcf/analysis/OE0538_projects_temp/DO-0008/data_temp/micromamba_test_rerun/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10")
data_hsc <- readRDS(snakemake@input[["data_hsc"]])
data <- data_hsc

print(data)

# prior markers determined
#markers_cons_hsc_prev <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/01_sign/signature_list_hsc")
#markers_cons_hsc <- readRDS("/omics/odcf/analysis/OE0538_projects_temp/DO-0008/data_temp/micromamba_test_rerun/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own/01_gens/geneset_list_hsc")
markers_cons_hsc <- readRDS(snakemake@input[["markers_cons_hsc"]])


metadata <- data@colData@listData %>% 
  as.data.frame() %>%
  select(-Sample)
rownames(metadata) <- paste(metadata$Barcode, metadata$Object_ID, sep = "_")

cleaner_data <- CreateSeuratObject(counts = data@assays@data@listData[["logcounts"]],
                                   meta.data = metadata)


# Identify markers that are specific to each cell type in each species
species <- unique(metadata$Species_ID)
cell_types <- unique(metadata$celltypes) %>% as.vector()
age <- unique(metadata$Age)

markers_list <- list()
markers_expression <- list()
for(y in age){
  for (s in species) {
    species_subset <- subset(cleaner_data,
                             subset = (Species_ID == s & Age == y))
    for (ct in cell_types) {
      # Get the cells for this cell type and species
      cells <- metadata %>%
        filter(Species_ID == s &
                 celltypes == ct &
                 Age == y) %>%
        row.names() %>% unique()
      if(length(cells) < 10){
        next()
      }
      # Identify markers that are specific to this cell type in this species
      markers <- FindMarkers(object = species_subset,
                             ident.1 = cells,
                             only.pos = TRUE,
                             min.pct = 0.1) # Require genes to be expressed in at least 10% of cells in each group
      print("done")
      # Add these markers to the list of markers for this cell type
      markers_list[[paste(s, ct, y, sep = "_")]] <- markers
      # identify what percentage of each animals' cells per type express identified markers
      animals <- grep(s,
                      unique(cleaner_data@meta.data[["Object_ID"]]),
                      value = T)
      expression <- matrix(data = NA,
                           nrow = nrow(markers),
                           ncol = length(animals))
      rownames(expression) <- rownames(markers)
      colnames(expression) <- animals
      for(a in animals){
        subset_animal <- subset(cleaner_data,
                                subset = Object_ID == as.character(a))
        expression[,a] <- suppressWarnings(
          DotPlot(subset_animal,
                  features = rownames(markers))$data[, "pct.exp"]
        )
      }

      markers_expression[[paste(s, ct, y, sep = "_")]] <- expression
    }
  }
}


markers_conservation_hsc <- list()
for(y in age){
  for(ct in cell_types){
    lists <- grep(paste(ct, y,  sep = "_"),
                  names(markers_list),
                  value = T)
    genes <- matrix(data = NA,
                    nrow = length(cleaner_data@assays[["RNA"]]@counts@Dimnames[[1]]),
                    ncol = length(species))
    rownames(genes) <- cleaner_data@assays[["RNA"]]@counts@Dimnames[[1]]
    colnames(genes) <- species
    for(l in lists){
      # remove gene if any of the individuals have 0 expression in the cell type
      expressed <- markers_expression[[l]][apply(markers_expression[[l]], 1,
                                                 function(x) !any(x == 0)), ]
      expression <- expressed %>% rowMeans() %>% as.data.frame()
      for(g in rownames(expression)){
        genes[g,substr(l,1,4)] <- expression[g,1]
      }
    }
    genes <- genes[rowSums(genes,
                           na.rm = T) > 0,]
    markers_conservation_hsc[[paste(ct, y,  sep = "_")]] <- genes
  }
}


#save(markers_conservation_hsc,
#     file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_hsc.RData")
save(markers_conservation_hsc,
     file = snakemake@output[["output_hsc"]])

# ## PLOT ##
# 
# load("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_hsc.RData")
# conserved_genes <- data.frame(celltype = rep(names(markers_cons_hsc)),
#                               young = 0,
#                               both = 0,
#                               old = 0,
#                               n_markers = 0,
#                               young_in_markers = 0,
#                               both_in_markers = 0,
#                               old_in_markers = 0)
# for(ct in names(markers_cons_hsc)){
#   young <- markers_conservation_hsc[[grep(paste(ct, "young", sep = "_"),
#                                       names(markers_conservation_hsc),
#                                       value = T)]] %>% as.data.frame()
#   young$gene <- rownames(young)
#   young <- young[complete.cases(young),]
# 
#   old <- markers_conservation_hsc[[grep(paste(ct, "old", sep = "_"),
#                                     names(markers_conservation_hsc),
#                                     value = T)]] %>% as.data.frame()
#   old$gene <- rownames(old)
#   old <- old[complete.cases(old),]
# 
#   both <- young$gene[young$gene %in% old$gene]
# 
#   mark <- markers_cons_hsc[[grep(ct,
#                                  names(markers_cons_hsc),
#                                  value = T)]][["conserved_df"]]
#   mark <- mark[complete.cases(mark),]
# 
#   conserved_genes[which(conserved_genes$celltype == ct),1] <- ct
# 
#   conserved_genes[which(conserved_genes$celltype == ct),2:8] <-
#     c(nrow(young),
#       length(both),
#       nrow(old),
#       nrow(mark),
#       sum(young$gene %in% mark$gene),
#       sum(both %in% mark$gene),
#       sum(old$gene %in% mark$gene)) %>% as.numeric()
# }


# plot_hsc <- data.frame(celltype = conserved_genes$celltype,
#                        both_in_mark = conserved_genes$both_in_markers,
#                        old_in_mark = (conserved_genes$both_in_markers/2) + (conserved_genes$old_in_markers - conserved_genes$both_in_markers),
#                        young_in_mark = -(conserved_genes$both_in_markers/2) - (conserved_genes$young_in_markers - conserved_genes$both_in_markers)) %>%
#   mutate(only_yng = young_in_mark - (conserved_genes$young - conserved_genes$young_in_markers),
#          only_old = old_in_mark + (conserved_genes$old - conserved_genes$old_in_markers),
#          both_in_mark_high = both_in_mark/2,
#          both_in_mark_low = -both_in_mark_high) %>%
#   select(-both_in_mark)
# 
# plot_hsc <- data.frame(celltype = plot_hsc$celltype,
#                        xmin = c(plot_hsc$only_yng, plot_hsc$young_in_mark,
#                                 plot_hsc$both_in_mark_low, plot_hsc$both_in_mark_high,
#                                 plot_hsc$old_in_mark),
#                        xmax = c(plot_hsc$young_in_mark,
#                                 plot_hsc$both_in_mark_low, plot_hsc$both_in_mark_high,
#                                 plot_hsc$old_in_mark, plot_hsc$only_old),
#                        fill = rep(c("red", "black", "gold", "black", "red"),
#                                   each = length(unique(plot_hsc$celltype))))
# 
# ggplot(plot_hsc,
#        aes(xmin = xmin,
#            xmax = xmax,
#            fill = fill)) +
#   facet_wrap(~celltype,
#              ncol = 1,
#              strip.position="left") +
#   geom_rect(aes(ymin = 0,
#                 ymax = 1)) +
#   theme_classic() +
#   theme(legend.position = "none",
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.background = element_blank()) +
#   scale_fill_manual(values = c("gold3", "gold", "grey40"))
# ggsave("02_age_figs/02_age_hsc.pdf",
#        height = length(names(markers_cons_hsc)),
#        width = 4)

# plot_hsc <- data.frame(celltype = conserved_genes$celltype,
#                        both_in_mark = conserved_genes$both_in_markers,
#                        old_and_mark = conserved_genes$old_in_markers - conserved_genes$both_in_markers,
#                        young_and_mark = conserved_genes$young_in_markers - conserved_genes$both_in_markers,
#                        only_yng = conserved_genes$young - conserved_genes$young_in_markers,
#                        only_old = conserved_genes$old - conserved_genes$old_in_markers) %>%
#   gather("category", "count", -celltype)
# 
# ggplot(plot_hsc,
#        aes(x = celltype,
#            y = count,
#            fill = factor(category,
#                          levels = c("only_yng",
#                                     "young_and_mark",
#                                     "both_in_mark",
#                                     "old_and_mark",
#                                     "only_old")))) +
#   geom_col(position = "fill") +
#   theme_classic() +
#   labs(x = NULL, y = "% marker genes identified", fill = NULL) +
#   theme(axis.text.x = element_text(angle = 90, h = 1))
# ggsave("02_age_figs/02_age_hsc_stack.pdf",
#        width = length(names(markers_cons_hsc)),
#        height = 4)
# 
# plot_hsc2 <- gather(conserved_genes[,c("celltype", "n_markers", "young_in_markers", "old_in_markers")], 
#                     "group", "count", -celltype)
# ggplot(plot_hsc2,
#        aes(x = celltype,
#            y = count,
#            fill = group)) +
#   geom_col(position = "dodge") +
#   theme_classic() +
#   labs(x = NULL, y = "# marker genes identified", fill = NULL) +
#   theme(axis.text.x = element_text(angle = 90, h = 1))
# ggsave("02_age_figs/02_age_hsc_totals.pdf",
#        width = length(names(markers_cons_hsc)),
#        height = 4)


###### Stromal

#data <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_str-10")
#markers_cons_str <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/01_sign/signature_list_str")

data_str <- readRDS(snakemake@input[["data_str"]])
data <- data_str

print(data)

markers_cons_str <- readRDS(snakemake@input[["markers_cons_str"]])

metadata <- data@colData@listData %>% 
  as.data.frame() %>%
  select(-Sample)
rownames(metadata) <- paste(metadata$Barcode, metadata$Object_ID, sep = "_")

cleaner_data <- CreateSeuratObject(counts = data@assays@data@listData[["logcounts"]],
                                   meta.data = metadata)


# Identify markers that are specific to each cell type in each species
species <- unique(metadata$Species_ID)
cell_types <- unique(metadata$celltypes) %>% as.vector()
age <- unique(metadata$Age)

markers_list <- list()
markers_expression <- list()
for(y in age){
  for (s in species) {
    species_subset <- subset(cleaner_data,
                             subset = (Species_ID == s & Age == y))
    for (ct in cell_types) {
      # Get the cells for this cell type and species
      cells <- metadata %>%
        filter(Species_ID == s &
                 celltypes == ct &
                 Age == y) %>%
        row.names() %>% unique()
      if(length(cells) < 10){
        next()
      }
      # Identify markers that are specific to this cell type in this species
      markers <- FindMarkers(object = species_subset,
                             ident.1 = cells,
                             only.pos = TRUE,
                             min.pct = 0.1) # Require genes to be expressed in at least 10% of cells in each group
      # Add these markers to the list of markers for this cell type
      markers_list[[paste(s, ct, y, sep = "_")]] <- markers
      # identify what percentage of each animals' cells per type express identified markers
      animals <- grep(s,
                      unique(cleaner_data@meta.data[["Object_ID"]]),
                      value = T)
      expression <- matrix(data = NA,
                           nrow = nrow(markers),
                           ncol = length(animals))
      rownames(expression) <- rownames(markers)
      colnames(expression) <- animals
      for(a in animals){
        subset_animal <- subset(cleaner_data,
                                subset = Object_ID == as.character(a))
        expression[,a] <- suppressWarnings(
          DotPlot(subset_animal,
                  features = rownames(markers))$data[, "pct.exp"]
        )
      }

      markers_expression[[paste(s, ct, y, sep = "_")]] <- expression
    }
  }
}


markers_conservation_str <- list()
for(y in age){
  for(ct in cell_types){
    lists <- grep(paste(ct, y,  sep = "_"),
                  names(markers_list),
                  value = T)
    genes <- matrix(data = NA,
                    nrow = length(cleaner_data@assays[["RNA"]]@counts@Dimnames[[1]]),
                    ncol = length(species))
    rownames(genes) <- cleaner_data@assays[["RNA"]]@counts@Dimnames[[1]]
    colnames(genes) <- species
    for(l in lists){
      # remove gene if any of the individuals have 0 expression in the cell type
      expressed <- markers_expression[[l]][apply(markers_expression[[l]], 1,
                                                 function(x) !any(x == 0)), ]
      expression <- expressed %>% rowMeans() %>% as.data.frame()
      for(g in rownames(expression)){
        genes[g,substr(l,1,4)] <- expression[g,1]
      }
    }
    genes <- genes[rowSums(genes,
                           na.rm = T) > 0,]
    markers_conservation_str[[paste(ct, y,  sep = "_")]] <- genes
  }
}

print(markers_conservation_str)

#save(markers_conservation_str,
#     file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_str.RData")
save(markers_conservation_str,
     file = snakemake@output[["output_str"]])

## PLOT ##

# load("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_str.RData")

# conserved_genes <- data.frame(celltype = rep(names(markers_cons_str)),
#                               young = 0,
#                               both = 0,
#                               old = 0,
#                               n_markers = 0,
#                               young_in_markers = 0,
#                               both_in_markers = 0,
#                               old_in_markers = 0)
# for(ct in names(markers_cons_str)){
#   young <- markers_conservation_str[[grep(paste(ct, "young", sep = "_"),
#                                       names(markers_conservation_str),
#                                       value = T)]] %>% as.data.frame()
#   young$gene <- rownames(young)
#   young <- young[complete.cases(young),]
# 
#   old <- markers_conservation_str[[grep(paste(ct, "old", sep = "_"),
#                                     names(markers_conservation_str),
#                                     value = T)]] %>% as.data.frame()
#   old$gene <- rownames(old)
#   old <- old[complete.cases(old),]
# 
#   both <- young$gene[young$gene %in% old$gene]
# 
#   mark <- markers_cons_str[[grep(ct,
#                              names(markers_cons_str),
#                              value = T)]][["conserved_df"]]
#   mark <- mark[complete.cases(mark),]
# 
#   conserved_genes[which(conserved_genes$celltype == ct),1] <- ct
# 
#   conserved_genes[which(conserved_genes$celltype == ct),2:8] <-
#     c(nrow(young),
#       length(both),
#       nrow(old),
#       nrow(mark),
#       sum(young$gene %in% mark$gene),
#       sum(both %in% mark$gene),
#       sum(old$gene %in% mark$gene)) %>% as.numeric()
# }


# plot_str <- data.frame(celltype = conserved_genes$celltype,
#                        both_in_mark = conserved_genes$both_in_markers,
#                        old_in_mark = (conserved_genes$both_in_markers/2) + (conserved_genes$old_in_markers - conserved_genes$both_in_markers),
#                        young_in_mark = -(conserved_genes$both_in_markers/2) - (conserved_genes$young_in_markers - conserved_genes$both_in_markers)) %>%
#   mutate(only_yng = young_in_mark - (conserved_genes$young - conserved_genes$young_in_markers),
#          only_old = old_in_mark + (conserved_genes$old - conserved_genes$old_in_markers),
#          both_in_mark_high = both_in_mark/2,
#          both_in_mark_low = -both_in_mark_high) %>%
#   select(-both_in_mark)
# plot_str <- data.frame(celltype = plot_str$celltype,
#                        xmin = c(plot_str$only_yng, plot_str$young_in_mark,
#                                 plot_str$both_in_mark_low, plot_str$both_in_mark_high,
#                                 plot_str$old_in_mark),
#                        xmax = c(plot_str$young_in_mark,
#                                 plot_str$both_in_mark_low, plot_str$both_in_mark_high,
#                                 plot_str$old_in_mark, plot_str$only_old),
#                        fill = rep(c("red", "black", "gold", "black", "red"),
#                                   each = length(unique(plot_str$celltype))))
# ggplot(plot_str,
#        aes(xmin = xmin,
#            xmax = xmax,
#            fill = fill)) +
#   facet_wrap(~celltype,
#              ncol = 1,
#              strip.position="left") +
#   geom_rect(aes(ymin = 0,
#                 ymax = 1)) +
#   theme_classic() +
#   theme(legend.position = "none",
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         strip.background = element_blank()) +
#   scale_fill_manual(values = c("gold3", "gold", "grey40"))
# ggsave("02_age_figs/02_age_str.pdf",
#        height = length(names(markers_cons_str)),
#        width = 4)

# plot_str <- data.frame(celltype = conserved_genes$celltype,
#                        both_in_mark = conserved_genes$both_in_markers,
#                        old_and_mark = conserved_genes$old_in_markers - conserved_genes$both_in_markers,
#                        young_and_mark = conserved_genes$young_in_markers - conserved_genes$both_in_markers,
#                        only_yng = conserved_genes$young - conserved_genes$young_in_markers,
#                        only_old = conserved_genes$old - conserved_genes$old_in_markers) %>%
#   gather("category", "count", -celltype)
# 
# ggplot(plot_str,
#        aes(x = celltype,
#            y = count,
#            fill = factor(category,
#                          levels = c("only_yng",
#                                     "young_and_mark",
#                                     "both_in_mark",
#                                     "old_and_mark",
#                                     "only_old")))) +
#   geom_col(position = "fill") +
#   theme_classic() +
#   labs(x = NULL, y = "% marker genes identified", fill = NULL) +
#   theme(axis.text.x = element_text(angle = 90, h = 1))
# ggsave("02_age_figs/02_age_str_stack.pdf",
#        width = length(names(markers_cons_str)),
#        height = 4)
# 
# plot_str2 <- gather(conserved_genes[,c("celltype", "n_markers", "young_in_markers", "old_in_markers")], 
#                     "group", "count", -celltype)
# ggplot(plot_str2,
#        aes(x = celltype,
#            y = count,
#            fill = group)) +
#   geom_col(position = "dodge") +
#   theme_classic() +
#   labs(x = NULL, y = "# marker genes identified", fill = NULL) +
#   theme(axis.text.x = element_text(angle = 90, h = 1))
# ggsave("02_age_figs/02_age_str_totals.pdf",
#        width = length(names(markers_cons_str)),
#        height = 4)
