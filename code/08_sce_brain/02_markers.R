library(Seurat)
library(dplyr)

data <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS")
data.updated <- UpdateSeuratObject(object = data)  # available data is v3 Seurat
DefaultAssay(data.updated) <- "RNA"

metadata <- data.updated@meta.data

# Identify markers that are specific to each cell type in each species
# ONLY LOOKING AT PRIMATES, WILL CROSS-COMPARE TO MOUSE
species <- c("human", "macaque", "marmoset") #unique(metadata$orig.ident)
cell_types <- unique(metadata$subclass_label)

markers_list <- list()
markers_expression <- list()
for (s in species) {
  species_subset <- subset(data.updated,
                           subset = orig.ident == s)
  for (ct in cell_types) {
    # Get the cells for this cell type and species
    cells <- metadata %>%
      filter(orig.ident == s &
               subclass_label == ct) %>%
      row.names() %>% unique()
    if(length(cells) ==  0){next()}
    # Identify markers that are specific to this cell type in this species
    markers <- FindMarkers(object = species_subset,
                           ident.1 = cells,
                           only.pos = TRUE,
                           min.pct = 0.1) # Require genes to be expressed in at least 10% of cells in each group
    # Add these markers to the list of markers for this cell type
    markers_list[[paste(s, ct, sep = "_")]] <- markers
    # identify what percentage of each animals' cells per type express identified markers
    animals <- unique(species_subset@meta.data[["donor"]])
    expression <- matrix(data = NA, 
                         nrow = nrow(markers), 
                         ncol = length(animals))
    rownames(expression) <- rownames(markers)
    colnames(expression) <- animals
    for(a in animals){
      subset_animal <- subset(data.updated,
                              subset = donor == as.character(a))
      expression[,a] <- suppressWarnings(
        DotPlot(subset_animal,
                features = rownames(markers),
                group.by = "orig.ident")$data[, "pct.exp"]
      )
    }
    
    markers_expression[[paste(s, ct, sep = "_")]] <- expression
  }
}

markers_conservation <- list()
for(ct in cell_types){
  lists <- grep(ct,
                names(markers_list),
                value = T)
  genes <- matrix(data = NA,
                  nrow = length(data.updated@assays[["RNA"]]@counts@Dimnames[[1]]),
                  ncol = length(species))
  rownames(genes) <- data.updated@assays[["RNA"]]@counts@Dimnames[[1]]
  colnames(genes) <- species
  for(l in lists){
    # remove gene if any of the individuals have 0 expression in the cell type
    expressed <- markers_expression[[l]][apply(markers_expression[[l]], 1, 
                                               function(x) !any(x == 0)), ]
    expression <- expressed %>% rowMeans() %>% as.data.frame()
    for(g in rownames(expression)){
      genes[g,strsplit(l, "_")[[1]][1]] <- expression[g,1]
    }
  }
  genes <- genes[rowSums(genes,
                         na.rm = T) > 0,]
  markers_conservation[[ct]] <- genes
}

saveRDS(markers_conservation,
        file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/08_sce_brain/02_marker_expression_primates.rds")

# since macaque is missing L2/3 IT cell type, first filter for non-all-NA columns
# then for non-NA rows to grab markers

core_markers <- list()
for(m in names(markers_conservation)){
  m_df <- markers_conservation[[m]]
  m_df <- m_df[, colSums(is.na(m_df)) < nrow(m_df)]
  m_df <- m_df[complete.cases(m_df),]
  
  core_markers[[m]] <- rownames(m_df)
}

saveRDS(core_markers,
        file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/08_sce_brain/02_marker_conserved_primates.rds")

