library(AUCell)
library(Seurat)
library(dplyr)
library(biomaRt)

setwd("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/04_rare_celltypes/mouse_lemur_Ezran2025")
mouse_lemur_harmony <- readRDS("01_mouse_lemur_harmony.rds")

hsc_sig <- read.delim("../../../../manuscript1/hsc_signature_table.csv", sep = ";", check.names = FALSE)
str_sig <- read.delim("../../../../manuscript1/str_signature_table.csv", sep = ";", check.names = FALSE)

sig_df <- bind_rows(hsc_sig, str_sig) %>% 
  dplyr::select(celltype, signature_gene)
gene_sets_mouse <- split(sig_df$signature_gene, sig_df$celltype)
#gene_sets_mouse <- lapply(gene_sets_mouse, unique)
names(gene_sets_mouse) <- make.names(names(gene_sets_mouse))

mmart <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")

mouse_symbols <- unique(unlist(gene_sets_mouse))
map <- getBM(attributes = c("external_gene_name",
                            "mmurinus_homolog_associated_gene_name",
                            "mmurinus_homolog_orthology_type"),
              filters    = "external_gene_name",
              values     = mouse_symbols,
              mart       = mmart)
colnames(map) <- c("mouse_symbol","lemur_symbol","orthology_type")
map <- map[map$lemur_symbol != "", ]
## actually erroneously removes some genes if do this filtering
# map <- map[order(map$orthology_type != "ortholog_one2one"), ]
# map <- map[!duplicated(map$mouse_symbol), ]

lemur_map <- setNames(map$lemur_symbol, map$mouse_symbol)
gene_sets_lemur <- lapply(gene_sets_mouse, function(g) unname(lemur_map[g]))
gene_sets_lemur <- lapply(gene_sets_lemur, function(v) v[!is.na(v)])

expr <- LayerData(mouse_lemur_harmony[["RNA"]], layer = "scale.data")
gene_sets_lemur <- lapply(gene_sets_lemur, function(g) intersect(g, rownames(expr)))

rankings <- AUCell_buildRankings(expr, 
                                 plotStats = FALSE, 
                                 verbose = FALSE)
cellsAUC  <- AUCell_calcAUC(gene_sets_lemur, 
                            rankings, 
                            aucMaxRank = ceiling(0.1 * nrow(rankings)))
auc_mat   <- t(getAUC(cellsAUC))

mouse_lemur_harmony <- AddMetaData(mouse_lemur_harmony, auc_mat)

saveRDS(mouse_lemur_harmony, 
  file = "02_mouse_lemur_harmony.rds")


