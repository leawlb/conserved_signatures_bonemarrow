library(AUCell)
library(Seurat)
library(dplyr)
library(biomaRt)

zebrafish_HSC <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/prepared/zeb_all_hspc")
hsc_sig <- read.delim("/omics/odcf/analysis/OE0538_projects/DO-0008/data/manuscript1/hsc_signature_table.csv", sep = ";", check.names = FALSE)

zebrafish_HSC <- zebrafish_HSC |>
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 10000) |>
  ScaleData() |>
  RunPCA(npcs = 50, verbose = FALSE) |>
  FindNeighbors(reduction = "pca", 
    dims = 1:50) |>
  RunUMAP(reduction = "pca", 
    dims = 1:50, 
    reduction.name = "umap_pca") |>
  FindClusters(resolution = 0.5)

DimPlot(zebrafish_HSC, 
  reduction = "umap_pca", 
  group.by = "cell_type", 
  raster = F) & NoAxes()


sig_df <- dplyr::select(hsc_sig, celltype, signature_gene)
gene_sets_mouse <- split(sig_df$signature_gene, sig_df$celltype)
names(gene_sets_mouse) <- make.names(names(gene_sets_mouse))

mmart <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")

mouse_symbols <- unique(unlist(gene_sets_mouse))
map <- getBM(attributes = c("external_gene_name",
                 "drerio_homolog_ensembl_gene",
                 "drerio_homolog_orthology_type"),
              filters    = "external_gene_name",
              values     = mouse_symbols,
              mart       = mmart)
colnames(map) <- c("mouse_symbol","zebrafish_symbol","orthology_type")
map <- map[map$zebrafish_symbol != "", ]

zebrafish_map <- setNames(map$zebrafish_symbol, map$mouse_symbol)
gene_sets_zebrafish <- lapply(gene_sets_mouse, function(g) unname(zebrafish_map[g]))
gene_sets_zebrafish <- lapply(gene_sets_zebrafish, function(v) v[!is.na(v)])


### stopped here, jumped below

expr <- GetAssayData(
  zebrafish_HSC,
  assay = "RNA",
  slot  = "scale.data"
)

gene_sets_zebrafish <- lapply(gene_sets_zebrafish, function(g) intersect(g, rownames(expr)))

rankings <- AUCell_buildRankings(expr, 
                                 plotStats = FALSE, 
                                 verbose = FALSE)
cellsAUC  <- AUCell_calcAUC(gene_sets_zebrafish, 
                            rankings, 
                            aucMaxRank = ceiling(0.1 * nrow(rankings)))
auc_mat   <- t(getAUC(cellsAUC))

zebrafish_HSC <- AddMetaData(zebrafish_HSC, auc_mat)

saveRDS(zebrafish_HSC, 
  file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/zebrafish_Athanasiadis2017/03_zebrafish.rds")







##### recieved data from original authors #####

meta <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/zebrafish_Athanasiadis2017/Tg_Tp_PhenoData.csv",
                  header = T, row.names = 1)
cells <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/zebrafish_Athanasiadis2017/All_Norm_Counts.csv",
                  header = T, row.names = 1)

rownames(cells) <- sub(",.*", "", rownames(cells))

ggplot(meta,
      aes(x = Component1,
          y = Component2,
          color = Type)) +
      geom_point() +
      theme_classic()


seu <- CreateSeuratObject(
  counts    = cells,      # we have no raw counts, so we use normalized
  meta.data = meta
)

expr <- GetAssayData(
  seu,
  assay = "RNA",
  slot  = "counts"
)
gene_sets_zebrafish <- lapply(gene_sets_zebrafish, function(g) intersect(g, rownames(expr)))

rankings <- AUCell_buildRankings(expr, 
                                 plotStats = FALSE, 
                                 verbose = FALSE)
cellsAUC  <- AUCell_calcAUC(gene_sets_zebrafish, 
                            rankings, 
                            aucMaxRank = ceiling(0.25 * nrow(rankings)))
auc_mat   <- t(getAUC(cellsAUC))

seu <- AddMetaData(seu, auc_mat)

saveRDS(seu, 
  file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/zebrafish_Athanasiadis2017/03_zebrafish_new.rds")

