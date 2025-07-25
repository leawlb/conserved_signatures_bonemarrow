library(Seurat)
library(dplyr)

########################################################
## identify stromal markers per species per cell type ##
########################################################

data_str <- readRDS(snakemake@input[["data_str"]])
data <- data_str

print(data)

metadata <- data@colData@listData %>% 
  as.data.frame() %>%
  select(-Sample)
rownames(metadata) <- paste(metadata$Barcode, metadata$Object_ID, sep = "_")

cleaner_data <- CreateSeuratObject(counts = data@assays@data@listData[["logcounts"]],
                                   meta.data = metadata)

# Identify markers that are specific to each cell type in each species
species <- unique(metadata$Species_ID)
cell_types <- unique(metadata$celltypes) %>% as.vector()

markers_list <- list()
markers_expression <- list()

for (s in species) {
  species_subset <- subset(cleaner_data,
                           subset = Species_ID == s)
  for (ct in cell_types) {
    # Get the cells for this cell type and species
    cells <- metadata %>%
      filter(Species_ID == s &
               celltypes == ct) %>%
      row.names() %>% unique()
    # Identify markers that are specific to this cell type in this species
    markers <- FindMarkers(object = species_subset,
                           ident.1 = cells,
                           only.pos = TRUE,
                           #logfc.threshold = 0.25, # default is 0.25
                           min.pct = 0.1) # Require genes to be expressed in at least 10% of cells in each group
    
    # Add these markers to the list of markers for this cell type
    markers_list[[paste(s, ct, sep = "_")]] <- markers
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
    
    markers_expression[[paste(s, ct, sep = "_")]] <- expression
  }
}

markers_conservation_str <- list()
for(ct in cell_types){
  lists <- grep(ct,
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
  markers_conservation_str[[ct]] <- genes
}

save(markers_conservation_str,
     file = snakemake@output[["output_str"]])

#########################
## same analysis, HSCs ##
#########################

print("starting HSCs")

data_hsc <- readRDS(snakemake@input[["data_hsc"]])
data <- data_hsc

print(data)

metadata <- data@colData@listData %>% 
  as.data.frame() %>%
  select(-Sample)
rownames(metadata) <- paste(metadata$Barcode, metadata$Object_ID, sep = "_")

cleaner_data <- CreateSeuratObject(counts = data@assays@data@listData[["logcounts"]],
                                   meta.data = metadata)

# Identify markers that are specific to each cell type in each species
species <- unique(metadata$Species_ID)
cell_types <- unique(metadata$celltypes) %>% as.vector()

markers_list <- list()
markers_expression <- list()
for (s in species) {
  species_subset <- subset(cleaner_data,
                           subset = Species_ID == s)
  for (ct in cell_types) {
    # Get the cells for this cell type and species
    cells <- metadata %>%
      filter(Species_ID == s &
               celltypes == ct) %>%
      row.names() %>% unique()
    # Identify markers that are specific to this cell type in this species
    markers <- FindMarkers(object = species_subset,
                           ident.1 = cells,
                           only.pos = TRUE,
                           #logfc.threshold = 0.25, # default is 0.25
                           min.pct = 0.1) # Require genes to be expressed in at least 10% of cells in each group
    # Add these markers to the list of markers for this cell type
    
    markers_list[[paste(s, ct, sep = "_")]] <- markers
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
    
    markers_expression[[paste(s, ct, sep = "_")]] <- expression
  }
}

markers_conservation_hsc <- list()
for(ct in cell_types){
  lists <- grep(ct,
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
  markers_conservation_hsc[[ct]] <- genes
}

save(markers_conservation_hsc,
     file = snakemake@output[["output_hsc"]])