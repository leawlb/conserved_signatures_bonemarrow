library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)

#base_path <- "/omics/odcf/analysis/OE0538_projects_temp/DO-0008/data_temp/micromamba_test_rerun/data"
#brain_path <- "/scRNAseq/main_analysis/sce_objects/08_sce_brain/"

#base_path_temp <- snakemake@params[["base_path"]] # try and see what happens if base_path is used --- I don't have access to OE0433
base_path <- snakemake@params[["base_path"]]
brain_path <- snakemake@params[["brain_path"]]

mouse_hs <- readRDS(paste0(base_path, brain_path, "04_recluster/mouse_reclust_hs.rds"))
mouse_sig <- readRDS(paste0(base_path, brain_path, "04_recluster/mouse_reclust_sig.rds"))

plot_hs <- mouse_hs@meta.data[,c("subclass_label", "seurat_clusters")] %>%
  table() %>% as.matrix()
plot_hs_quant <- plot_hs %>%
  as.data.frame() %>%
  group_by(seurat_clusters) %>%
  mutate(clust_prop = Freq/sum(Freq),
         clust_assigned = (max(clust_prop)==clust_prop)) %>%
  filter(subclass_label != "L6 IT Car3") %>%
  group_by(clust_assigned, subclass_label) %>%
  mutate(Count = sum(Freq)) %>%
  select(-seurat_clusters, -clust_prop, -Freq) %>%
  unique()
sum(filter(plot_hs_quant, clust_assigned == T)$Count)/sum(plot_hs_quant$Count)
# [1] 0.9060639

plot_sig <- mouse_sig@meta.data[,c("subclass_label", "seurat_clusters")] %>%
  table() %>% as.matrix()
plot_sig_quant <- plot_sig %>%
  as.data.frame() %>%
  group_by(seurat_clusters) %>%
  mutate(clust_prop = Freq/sum(Freq),
         clust_assigned = (max(clust_prop)==clust_prop)) %>%
  filter(subclass_label != "L6 IT Car3") %>%
  group_by(clust_assigned, subclass_label) %>%
  mutate(Count = sum(Freq)) %>%
  select(-seurat_clusters, -clust_prop, -Freq) %>%
  unique()
sum(filter(plot_sig_quant, clust_assigned == T)$Count)/sum(plot_sig_quant$Count)
# [1] 0.9172439

all_data <- rbind(
  mutate(plot_hs_quant, genes = "human markers"),
  mutate(plot_sig_quant, genes = "signature")
)

ggplot(all_data,
       aes(y = factor(genes, 
                      levels = c("signature", "human markers")),
           x = Count,
           fill = clust_assigned)) +
  facet_wrap(~factor(subclass_label, 
                     levels = c("L2/3 IT", "L5 IT", "L6 IT", 
                                "L5 ET", "L5/6 NP", "L6 CT", "L6b")), 
             ncol = 1,
             strip.position = "right") +
  geom_col(position = "stack") +
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0),
    strip.background = element_rect(color = "white")) +
  labs(y = NULL, 
       x = "Mouse cells", 
       fill = "Accurately\nAssigned") +
  scale_fill_manual(values = c("grey90", "blue"))

########################
## unused UMAP visual ##
########################

plot_hs_clust <- plot_hs %>%
  as.data.frame() %>%
  group_by(seurat_clusters) %>%
  mutate(clust_assigned = (max(Freq)==Freq)) %>%
  filter(clust_assigned == T) %>%
  select(-clust_assigned, -Freq)
test <- lapply(mouse_hs@meta.data$seurat_clusters, 
               function(x){
                 plot_hs_clust[which(plot_hs_clust$seurat_clusters == x), "subclass_label"]}) %>% 
  unlist()
mouse_hs@meta.data$correct <- mouse_hs@meta.data$subclass_label == test

plot_sig_clust <- plot_sig %>%
  as.data.frame() %>%
  filter(subclass_label != "L6 IT Car3") %>%
  group_by(seurat_clusters) %>%
  mutate(clust_assigned = (max(Freq)==Freq)) %>%
  filter(clust_assigned == T) %>%
  select(-clust_assigned, -Freq)
test <- lapply(mouse_sig@meta.data$seurat_clusters, 
               function(x){
                 plot_sig_clust[which(plot_sig_clust$seurat_clusters == x), "subclass_label"]}) %>% 
  unlist()
mouse_sig@meta.data$correct <- mouse_sig@meta.data$subclass_label == test

rownames <- rownames(mouse_sig@meta.data)
mouse_sig@meta.data <- merge(mouse_sig@meta.data, # correct.x is signature markers
                             mouse_hs@meta.data[,c("nCount_RNA", "nFeature_RNA", "sample_id", "correct")], # correct.y is human markers
                             by = c("nCount_RNA", "nFeature_RNA", "sample_id")) %>%
  mutate(correct = paste(correct.x, correct.y))
rownames(mouse_sig@meta.data) <- rownames

# > table(mouse_sig@meta.data$correct)
# 
# FALSE FALSE  FALSE TRUE  TRUE FALSE   TRUE TRUE 
# 388         306         398        7150 

DimPlot(mouse_sig,
        group.by = "correct") +
  theme_void() +
  labs(title = NULL)


