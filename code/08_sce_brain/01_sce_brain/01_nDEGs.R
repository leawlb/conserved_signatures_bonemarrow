library(tidyverse, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(scuttle, quietly = TRUE)
library(BiocGenerics, quietly = TRUE)
set.seed(37)

# Load and update Seurat object
#data <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/08_sce_brain/sample.combined_exc_4_species_integration.RDS")
data <- readRDS(snakemake@input[["data_input"]])

data.updated <- UpdateSeuratObject(object = data)  # available data is v3 Seurat
DefaultAssay(data.updated) <- "RNA"

# Convert Seurat object to SingleCellExperiment
data.sce <- as.SingleCellExperiment(data.updated)

# Aggregate counts across cells by donor and subclass_label
agg <- scuttle::aggregateAcrossCells(
  data.sce,
  use.assay.type = "counts",
  id = colData(data.sce)[, c("donor", "subclass_label")]
)

# Split aggregated data by subclass_label
agg_list <- list()
for (i in unique(agg$subclass_label)) {
  agg_list[[i]] <- agg[, agg$subclass_label == i]
  agg_list[[i]]$sample <- factor(agg_list[[i]]$donor)
  
  colData(agg_list[[i]]) <- colData(agg_list[[i]])[, c("orig.ident", 
                                                       "subclass_label", 
                                                       "donor")]
}

# Define combinations for contrasts
comb_list <- list(
  c("human", "macaque"),
  c("human", "marmoset"),
  c("macaque", "marmoset")
)

results_list <- list()
# Loop over each subclass label (dataset) in agg_list
for (subclass in names(agg_list)) {
  se <- agg_list[[subclass]]
  se$orig.ident <- as.factor(se$orig.ident)
  
  # Initialize list to store results for this dataset
  pairwise_results <- list()
  
  # Loop over each combination in comb_list
  for (comb in comb_list) {
    if (all(comb %in% levels(se$orig.ident))) {
      # Set up DESeqDataSet
      dds <- DESeqDataSet(se, design = ~ orig.ident)
      dds <- DESeq(dds)
      
      # Define the contrast
      contrast <- c("orig.ident", comb[1], comb[2])
      
      # Extract results
      res <- results(dds, contrast = contrast)
      res_filtered <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
      
      # Store results for this combination
      pairwise_results[[paste(comb, collapse = "_vs_")]] <- res_filtered
    } else {
      warning(paste("Combination", paste(comb, collapse = " vs "), "is not fully represented in the data."))
      pairwise_results[[paste(comb, collapse = "_vs_")]] <- NULL
    }
  }
  
  # Store results for this subclass label
  results_list[[subclass]] <- pairwise_results
}

# Set cutoff values
fc_cutoff <- 1.5
padj_cutoff <- 0.05

# make a list of genes shared between each pairwise comp, for each cell type
res_list_shared_comp <- lapply(results_list, function(pairwise_results){
  
  shared_list <- list()
  
  # Extract comparison names
  comparisons <- names(pairwise_results)
  
  for (comp_name in comparisons){
    res_df <- pairwise_results[[comp_name]]
    
    if (!is.null(res_df)){
      # Get genes with absolute log2 FC < cutoff and padj > cutoff
      shared_genes <- rownames(res_df)[
        abs(res_df$log2FoldChange) < fc_cutoff & 
          res_df$padj > padj_cutoff]
      
      # Export unique shared genes
      shared_genes <- unique(shared_genes)
      shared_list[[comp_name]] <- shared_genes
    } else {
      shared_list[[comp_name]] <- character(0)
    }
  }
  
  return(shared_list)
})


# make a list of genes that are shared between ALL pairwise comparisons
res_list_shared_all <- lapply(res_list_shared_comp, function(shared){
  
  shared_all <- BiocGenerics::Reduce(BiocGenerics::intersect, shared)
  
  return(shared_all)
})

#saveRDS(res_list_shared_all, 
#        file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/08_sce_brain/01_list_nDEGs_all.rds")

saveRDS(res_list_shared_all, snakemake@output[["data_output"]])
