#-------------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(scater, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------

fraction_curr <- snakemake@wildcards[["fraction"]]
print(fraction_curr)

signature_list <- base::readRDS(snakemake@input[["signature_list"]])
sce <- base::readRDS(snakemake@input[["sce_inp"]])

ensembl_hum <- base::readRDS(snakemake@input[["ensembl_hum"]])
ensembl_mus <- base::readRDS(snakemake@input[["ensembl_mus"]])
#ensembl_mus <- base::readRDS(snakemake@input[["ensembl_hgl"]])
#ensembl_mus <- base::readRDS(snakemake@input[["ensembl_rdn"]])

#-------------------------------------------------------------------------------

# get the unique list of all signature genes, conserved marker genes, or ndges 
# per fraction
# get the IDs from the SCE object rowData, because this is where
# the gene symbol came from originally
sign <- vector()
for(i in 1:length(signature_list)){
  sign <- c(sign, signature_list[[i]]$conserved_signature)
}
sign <- base::unique(sign)

mark <- vector()
for(i in 1:length(signature_list)){
  mark <- c(mark, signature_list[[i]]$conserved_markers)
}
mark <- base::unique(mark)

ndge <- vector()
for(i in 1:length(signature_list)){
  ndge <- c(ndge, signature_list[[i]]$ndges)
}
ndge <- base::unique(ndge)

subclustering_genes <- vector()
for(i in 1:length(signature_list)){
  subclustering_genes <- c(subclustering_genes, signature_list[[i]]$genes_subclustering)
}
subclustering_genes <- base::unique(subclustering_genes)

# remove subclustering genes from gene lists because they were manually picked
consm <- consm[!consm %in% subclustering_genes]
ndges <- ndges[!ndges %in% subclustering_genes]
signt <- signt[!signt %in% subclustering_genes]

# get the IDs
sign_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% sign]
mark_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% mark]
ndge_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% ndge]

length(sign_ids)
length(mark_ids)
length(ndge_ids)

mouse_ids <- base::unique(c(sign_ids, mark_ids, ndge_ids))

#-------------------------------------------------------------------------------
# TODO: change it to add hgl and dnr later
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
human_ids1 <- base::unique(ensembl_mus_sub$hsapiens_homolog_ensembl_gene)

# from human ensembl dataframe (via musmusculus homologue IDs)
human_ids2 <- base::unique(ensembl_hum$ensembl_gene_id[
  ensembl_hum$mmusculus_homolog_ensembl_gene %in% mouse_ids])

base::table(is.na(base::match(human_ids1, human_ids2)))
base::table(is.na(base::match(human_ids2, human_ids1)))

# use the human homologue IDs which were already filtered above
human_ids <- human_ids1
ensembl_hum_sub <- ensembl_hum %>%
  dplyr::filter(ensembl_gene_id %in% human_ids)

# remove duplicated genes that have multiple corresponding mouse IDs
ensembl_hum_sub <- ensembl_hum_sub %>% 
  dplyr::filter(!ensembl_gene_id %in% 
                  base::unique(ensembl_gene_id[
                    base::duplicated(ensembl_gene_id)])) 

# add mouse gene name
ensembl_hum_sub$mmusculus_gene_name <- ensembl_mus_sub$external_gene_name[
  base::match(ensembl_hum_sub$ensembl_gene_id, 
              ensembl_mus_sub$hsapiens_homolog_ensembl_gene)]

ensembl_hum_sub$mmusculus_gene_name_sce <- rowData(sce)$Symbol[
  base::match(ensembl_hum_sub$mmusculus_homolog_ensembl_gene, 
              rowData(sce)$ID)]

base::table(base::duplicated(ensembl_hum_sub$ensembl_gene_id))
base::table(base::duplicated(ensembl_hum_sub$external_gene_name))
base::table(base::duplicated(ensembl_hum_sub$mmusculus_homolog_ensembl_gene))
base::table(base::duplicated(ensembl_hum_sub$mmusculus_gene_name))
base::table(base::duplicated(ensembl_hum_sub$mmusculus_gene_name_sce))

table(is.na(base::match(ensembl_hum_sub$mmusculus_external_gene_name,
                        ensembl_hum_sub$mmusculus_gene_name_sce)))

head(ensembl_hum_sub)

#-------------------------------------------------------------------------------
# export different sets

ensembl_emfs <- ensembl_hum_sub %>%
  dplyr::filter(mmusculus_homolog_ensembl_gene %in% emfs_ids)
ensembl_mark <- ensembl_hum_sub %>%
  dplyr::filter(mmusculus_homolog_ensembl_gene %in% mark_ids)
ensembl_ndge <- ensembl_hum_sub %>%
  dplyr::filter(mmusculus_homolog_ensembl_gene %in% ndge_ids)

base::saveRDS(ensembl_emfs, snakemake@output[["ensembl_sign"]])
base::saveRDS(ensembl_mark, snakemake@output[["ensembl_mark"]])
base::saveRDS(ensembl_ndge, snakemake@output[["ensembl_ndge"]])
