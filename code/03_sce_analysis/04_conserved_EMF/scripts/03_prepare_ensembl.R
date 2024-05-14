#-------------------------------------------------------------------------------

library(tidyverse)
library(scater)
set.seed(37)

#-------------------------------------------------------------------------------

fraction_curr <- snakemake@wildcards[["fraction"]]
print(fraction_curr)

# fraction_curr <- "hsc"
#sce <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_str-10")
#emf_list <- readRDS(paste0("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_conserved_EMF/01_emfs/emf_list_", fraction_curr))
# sce <- readRDS(paste0("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/02_sce_anno/10_anns/sce_", fraction_curr, "-10"))
#ensembl_hum <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/ensembl_hum")
#ensembl_mus <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/ensembl_mus")

emf_list <- readRDS(snakemake@input[["emf_list_inp"]])
sce <- readRDS(snakemake@input[["sce_inp"]])

ensembl_hum <- readRDS(snakemake@input[["ensembl_hum"]])
ensembl_mus <- readRDS(snakemake@input[["ensembl_mus"]])

#-------------------------------------------------------------------------------

# get the unique list of all emfs, conserved marker genes, or ndges per fraction

emfs <- vector()
for(i in 1:length(emf_list)){
  emfs <- c(emfs, emf_list[[i]]$cons_emf)
}
emfs <- unique(emfs)
emfs_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% emfs]

mark <- vector()
for(i in 1:length(emf_list)){
  mark <- c(mark, emf_list[[i]]$df_cons$gene)
}
mark <- unique(mark)
mark_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% mark]

ndge <- vector()
for(i in 1:length(emf_list)){
  ndge <- c(ndge, emf_list[[i]]$ndges)
}
ndge <- unique(ndge)
ndge_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% ndge]

length(emfs_ids)
length(mark_ids)
length(ndge_ids)

mouse_ids <- unique(c(emfs_ids, mark_ids, ndge_ids))

#-------------------------------------------------------------------------------

# get the corresponding IDs from the SCE object and subset ensembl_mus
ensembl_mus_sub <- ensembl_mus %>% 
  dplyr::filter(ensembl_gene_id %in% mouse_ids)

# get duplicated mouse ensembl IDs and show all
dup_mouse_ids <- ensembl_mus_sub %>% 
  dplyr::filter(duplicated(ensembl_gene_id)) %>%
  dplyr::pull(ensembl_gene_id)

ensembl_mus_sub %>% 
  dplyr::filter(ensembl_gene_id %in% dup_mouse_ids)

# remove all duplicated gene IDs 
# = mouse IDs with multiple human homologue IDs
# = cannot uniquely identify one human gene
ensembl_mus_sub <- ensembl_mus_sub %>%
  dplyr::filter(!ensembl_gene_id %in% dup_mouse_ids)

head(ensembl_mus_sub)

#-------------------------------------------------------------------------------

# check which human IDS correspond to mouse IDs
# from mouse ensembl dataframe (via human homologue IDs)
human_ids1 <- unique(ensembl_mus_sub$hsapiens_homolog_ensembl_gene)

# from human ensembl dataframe (via musmusculus homologue IDs)
human_ids2 <- unique(ensembl_hum$ensembl_gene_id[
  ensembl_hum$mmusculus_homolog_ensembl_gene %in% mouse_ids])

table(is.na(match(human_ids1, human_ids2)))
table(is.na(match(human_ids2, human_ids1)))

# use the human homologue IDs which were already filtered above
human_ids <- human_ids1
ensembl_hum_sub <- ensembl_hum %>%
  dplyr::filter(ensembl_gene_id %in% human_ids)

# remove duplicated genes that have multiple corresponding mouse IDs
ensembl_hum_sub <- ensembl_hum_sub %>% 
  dplyr::filter(!ensembl_gene_id %in% 
           unique(ensembl_gene_id[
             duplicated(ensembl_gene_id)])) 

# add mouse gene name
ensembl_hum_sub$mmusculus_gene_name <- ensembl_mus_sub$external_gene_name[
  match(ensembl_hum_sub$ensembl_gene_id, 
        ensembl_mus_sub$hsapiens_homolog_ensembl_gene)]

ensembl_hum_sub$mmusculus_gene_name_sce <- rowData(sce)$Symbol[
  match(ensembl_hum_sub$mmusculus_homolog_ensembl_gene, 
        rowData(sce)$ID)]

table(duplicated(ensembl_hum_sub$ensembl_gene_id))
table(duplicated(ensembl_hum_sub$external_gene_name))
table(duplicated(ensembl_hum_sub$mmusculus_homolog_ensembl_gene))
table(duplicated(ensembl_hum_sub$mmusculus_gene_name))
table(duplicated(ensembl_hum_sub$mmusculus_gene_name_sce))

table(is.na(match(ensembl_hum_sub$mmusculus_external_gene_name,
                  ensembl_hum_sub$mmusculus_gene_name_sce)))

head(ensembl_hum_sub)

#-------------------------------------------------------------------------------
# export different sets

ensembl_emfs <- ensembl_hum_sub %>%
  filter(mmusculus_homolog_ensembl_gene %in% emfs_ids)
ensembl_mark <- ensembl_hum_sub %>%
  filter(mmusculus_homolog_ensembl_gene %in% mark_ids)
ensembl_ndge <- ensembl_hum_sub %>%
  filter(mmusculus_homolog_ensembl_gene %in% ndge_ids)

saveRDS(ensembl_emfs, snakemake@output[["ensembl_emfs"]])
saveRDS(ensembl_mark, snakemake@output[["ensembl_mark"]])
saveRDS(ensembl_ndge, snakemake@output[["ensembl_ndge"]])
