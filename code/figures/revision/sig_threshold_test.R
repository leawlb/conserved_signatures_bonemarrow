library(tidyverse)
library(Seurat)

load("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_hsc.RData")

hsc_nDGEs <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/02_DESeq2_crossspecies/05_nres/PC_0.05_FC_1.5/res_hsc_celltype_dfs")

hsc_sig <- read.delim("/omics/odcf/analysis/OE0538_projects/DO-0008/data/manuscript1/hsc_signature_table.csv", sep = ";")

setwd("/omics/groups/OE0433/internal/veronica/github/Interspecies_BM_phd/")

################################
##### modify pval of nDEGs #####
################################

p_try <- c(0.000001, 0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.1)

p_chr <- as.character(p_try)
signature_genes_by_p <- setNames(vector("list", length(p_chr)), p_chr)

res_list <- list()
i <- 1L

for (ct in names(hsc_nDGEs)) {
  expected <- rownames(markers_conservation_hsc[[ct]][
    complete.cases(markers_conservation_hsc[[ct]]), 
  ])

  for (p in p_try) {
    genes_p <- hsc_nDGEs[[ct]] %>%
      dplyr::filter(padj < p) %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(dplyr::n() == 12) %>%   # 12 pairwise contrasts
      dplyr::pull(gene) %>%
      unique()

    overlap_genes <- intersect(expected, genes_p)

    res_list[[i]] <- tibble::tibble(
      celltype = ct,
      padj = p,
      n_nDEGs = length(genes_p),
      n_signature_genes = length(overlap_genes)
    )
    i <- i + 1L

    # accumulate unique signature genes per p across cell types
    key <- as.character(p)
    signature_genes_by_p[[key]] <- union(signature_genes_by_p[[key]], overlap_genes)
  }
}

nDEG_df <- dplyr::bind_rows(res_list)

saveRDS(nDEG_df, "nDEG_df.rds")
saveRDS(signature_genes_by_p, "signature_genes_by_p.rds")

scale_factor <- max(nDEG_df$n_nDEGs, na.rm = TRUE) /
                max(nDEG_df$n_signature_genes, na.rm = TRUE)
pdf(file = "code/figures/revision/sig_threshold_DEGs.pdf",
    height = 7, width = 7)
ggplot(nDEG_df, 
  aes(x = -log(padj), 
      color = celltype)) +
  geom_vline(xintercept = -log(0.05), 
              color = "grey80") +
  geom_line(aes(y = final_count),   # left-axis data
                linetype = "solid") +
  geom_line(aes(y = total_count / scale_factor),  # right-axis data, rescaled
                linetype = "dashed") +
  scale_y_continuous(
    name = "Total signature genes (solid)",
    sec.axis = sec_axis(~ . * scale_factor,
                        name = "Total nDEGs (dotted)")) +
  labs(color = NULL, x = "-log(p.adj) nDEG analysis") +
  theme_classic()
dev.off()




################################
#### modify l2FC of markers ####
################################

data <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_hsc-10")

t_try <- c(0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1)

metadata <- data@colData@listData %>%
  as.data.frame() %>%
  select(-Sample)

rownames(metadata) <- paste(metadata$Barcode, metadata$Object_ID, sep = "_")

cleaner_data <- CreateSeuratObject(
  counts    = data@assays@data@listData[["logcounts"]],
  meta.data = metadata
)

species    <- unique(metadata$Species_ID)
cell_types <- unique(metadata$celltypes)

## nDGE signature genes per cell type (padj < 0.05 and present in all 12 contrasts)
ndge_genes_by_ct <- setNames(
  lapply(cell_types, function(ct) {
    x <- hsc_nDGEs[[ct]]
    if (is.null(x)) return(character(0))

    x %>%
      filter(padj < 0.05) %>%
      count(gene) %>%
      filter(n == 12) %>%
      pull(gene) %>%
      unique()
  }),
  cell_types
)

## containers
t_chr <- as.character(t_try)
markers_all_thresholds <- vector("list", length(t_chr)); names(markers_all_thresholds) <- t_chr
common_marker_stats    <- vector("list", length(t_chr)); names(common_marker_stats)    <- t_chr
signature_genes_by_t   <- vector("list", length(t_chr)); names(signature_genes_by_t)   <- t_chr

for (t in t_try) {
  message("logfc.threshold = ", t)
  t_name <- as.character(t)

  ## markers for this threshold: species × cell type
  markers_list <- list()

  for (s in species) {
    species_subset <- subset(cleaner_data, subset = Species_ID == s)

    for (ct in cell_types) {
      cells <- metadata %>%
        filter(Species_ID == s, celltypes == ct) %>%
        rownames()

      cells_in_obj <- intersect(colnames(species_subset), cells)
      if (length(cells_in_obj) < 3L) next

      markers <- FindMarkers(
        object          = species_subset,
        ident.1         = cells_in_obj,
        only.pos        = TRUE,
        logfc.threshold = t,
        min.pct         = 0.1
      )

      if (nrow(markers) == 0L) next

      markers_list[[paste(s, ct, sep = "_")]] <- markers
    }
  }

  markers_all_thresholds[[t_name]] <- markers_list

  ## stats + overlaps for this threshold
  stats_t <- bind_rows(lapply(cell_types, function(ct) {
    # markers per species for this cell type
    genes_by_species <- lapply(species, function(s) {
      key <- paste(s, ct, sep = "_")
      if (is.null(markers_list[[key]])) character(0) else rownames(markers_list[[key]])
    })

    # require genes in all species
    if (any(lengths(genes_by_species) == 0L)) {
      common_genes <- character(0)
    } else {
      common_genes <- Reduce(intersect, genes_by_species)
    }

    overlap <- intersect(common_genes, ndge_genes_by_ct[[ct]])

    tibble(
      celltype                     = ct,
      logfc_threshold              = t,
      n_common_markers_all_species = length(common_genes),
      n_signature_genes            = length(overlap),
      signature_genes              = list(overlap)
    )
  }))

  common_marker_stats[[t_name]]  <- select(stats_t, -signature_genes)
  signature_genes_by_t[[t_name]] <- unique(unlist(stats_t$signature_genes))
}

## final summary table across all thresholds
common_marker_stats_df <- bind_rows(common_marker_stats)

pdf(file = "code/figures/revision/sig_threshold_markers.pdf",
    height = 7, width = 7)
ggplot(common_marker_stats_df, 
  aes(x = logfc_threshold, 
      color = celltype)) +
  geom_vline(xintercept = 0.25, 
              color = "grey80") +
  geom_line(aes(y = n_signature_genes),
                linetype = "solid") +
  geom_line(aes(y = n_common_markers_all_species),
                linetype = "dashed") +
  scale_y_log10(
    name = "Total signature genes (solid),\nTotal conserved markers (dotted)") +
  labs(color = NULL, x = "marker log2FC") +
  theme_classic()
dev.off()


saveRDS(markers_all_thresholds, "markers_all_thresholds.rds")
saveRDS(common_marker_stats_df, "common_marker_stats_df.rds")
saveRDS(signature_genes_by_t, "signature_genes_by_t.rds")
