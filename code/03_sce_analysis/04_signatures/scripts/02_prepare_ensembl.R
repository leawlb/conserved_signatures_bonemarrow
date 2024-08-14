#-------------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
set.seed(37)

#-------------------------------------------------------------------------------

fraction_curr <- snakemake@wildcards[["fraction"]]
print(fraction_curr)

signature_list <- base::readRDS(snakemake@input[["signature_list"]])
sce <- base::readRDS(snakemake@input[["sce_inp"]])

ensembl_mus <- base::readRDS(snakemake@input[["ensembl_mus"]])
ensembl_hum <- base::readRDS(snakemake@input[["ensembl_hum"]])
ensembl_zeb <- base::readRDS(snakemake@input[["ensembl_zeb"]])
ensembl_nmr <- base::readRDS(snakemake@input[["ensembl_nmr"]])

print(head(ensembl_mus))
print(head(ensembl_hum))
print(head(ensembl_zeb))
print(head(ensembl_nmr))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# prepare gene sets

# get the unique list of all signature genes, conserved marker genes, or ndges 
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

mmms <- vector()
for(i in 1:length(signature_list)){
  mmms <- c(mmms, signature_list[[i]]$conserved_df$gene[
    which(!is.na(signature_list[[i]]$conserved_df$mmus))])
}
mmms <- base::unique(mmms)

# subclustering_genes <- vector()
# for(i in 1:length(signature_list)){
#   subclustering_genes <- c(subclustering_genes,
#                            signature_list[[i]]$genes_subclustering)
# }
# subclustering_genes <- base::unique(subclustering_genes)

# remove subclustering genes from gene lists because they were manually picked
# sign <- sign[!sign %in% subclustering_genes]
# mark <- mark[!mark %in% subclustering_genes]
# ndge <- ndge[!ndge %in% subclustering_genes]
# mmms <- mmms[!mmms %in% subclustering_genes]
# for other species it is alright to keep the genes even if they were used for
# subclustering of our dataset - they were still identified as part of a set

#-------------------------------------------------------------------------------
# get the IDs from the SCE object rowData, this is where
# the gene symbol came from originally
sign_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% sign]
mark_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% mark]
ndge_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% ndge]
mmms_ids <- rowData(sce)$ID[rowData(sce)$Symbol %in% mmms]

length(sign_ids)
length(mark_ids)
length(ndge_ids)
length(mmms_ids)

mouse_ids <- base::unique(c(sign_ids, mark_ids, ndge_ids, mmms_ids))
print(length(mouse_ids))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# prepare ensembl dataframe for all species

ensembl_df <- ensembl_mus

# unify colnames
colnames(ensembl_df) <- c("ENSMUS_ID", "MMUS_SYMBOL")
colnames(ensembl_hum) <- c("ENSG_ID", "HSAPIENS_SYMBOL", "ENSMUS_ID")
colnames(ensembl_zeb) <- c("ENSDARG_ID", "DRERIO_SYMBOL", "ENSMUS_ID")
colnames(ensembl_nmr) <- c("ENSHGL_ID", "HGLABER_SYMBOL", "ENSMUS_ID")

# shorten them so the ensembl_df doesn't become too large
ensembl_hum <- ensembl_hum[ensembl_hum$ENSMUS_ID %in% mouse_ids,]
ensembl_zeb <- ensembl_zeb[ensembl_zeb$ENSMUS_ID %in% mouse_ids,]
ensembl_nmr <- ensembl_nmr[ensembl_nmr$ENSMUS_ID %in% mouse_ids,]

# join all dataframes based on mouse IDs
ensembl_df <- ensembl_df %>%  
  dplyr::left_join(ensembl_hum,
                   join_by(ENSMUS_ID),
                   relationship = "many-to-many") |>
  dplyr::left_join(ensembl_zeb,
                   join_by(ENSMUS_ID),
                   relationship = "many-to-many") |>
  dplyr::left_join(ensembl_nmr,
                   join_by(ENSMUS_ID),
                   relationship = "many-to-many") |>
  dplyr::filter(ENSMUS_ID %in% mouse_ids)

print(base::table(base::duplicated(ensembl_df$ENSMUS_ID)))
print(head(ensembl_df))

# get duplicated mouse ensembl IDs and show all
# these are genes with multiple species homologue IDs
dup_mouse_ids <- ensembl_df %>% 
  dplyr::filter(base::duplicated(ENSMUS_ID)) %>%
  dplyr::pull(ENSMUS_ID) %>%
  base::unique()
#print(dup_mouse_ids)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# subset ensembl_df to required gene sets

ensembl_df_sign <- ensembl_df %>%
  dplyr::filter(ENSMUS_ID %in% sign_ids)
ensembl_df_mark <- ensembl_df %>%
  dplyr::filter(ENSMUS_ID %in% mark_ids)
ensembl_df_ndge <- ensembl_df %>%
  dplyr::filter(ENSMUS_ID %in% ndge_ids)
ensembl_df_mmms <- ensembl_df %>%
  dplyr::filter(ENSMUS_ID %in% mmms_ids)

base::saveRDS(ensembl_df_sign, snakemake@output[["ensembl_sign"]])
base::saveRDS(ensembl_df_mark, snakemake@output[["ensembl_mark"]])
base::saveRDS(ensembl_df_ndge, snakemake@output[["ensembl_ndge"]])
base::saveRDS(ensembl_df_mmms, snakemake@output[["ensembl_mmms"]])
