
mcar <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_mrge/sce_mcar-07")
mcas <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_ctyp/sce_mcas-07")
mmus <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_ctyp/sce_mmus-07")
mspr <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/07_ctyp/sce_mspr-07")

sce_dol <- readRDS(file =  "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/ref_dolgalev_sce")
sce_dah <- readRDS(file =  "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/ref_dahlin_sce")
sce_bac <- readRDS(file =  "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/ref_baccin_sce")


as.vector(table(mmus$Identity_ref_fraction[mmus$Age_ID == "yng"]))
names(table(mmus$Identity_ref_fraction[mmus$Age_ID == "yng"]))

# Dolgalev
celltype_df_dol <- data.frame(
  celltype = names(table(sce_dol$identity_ref))
)

celltype_df_dol$mmus_yng <- as.vector(
  table(mmus$Identity_ref_fraction[mmus$Age_ID == "yng"]))[
    match(celltype_df_dol$celltype, names(table(
      mmus$Identity_ref_fraction[mmus$Age_ID == "yng"])))]

celltype_df_dol$mmus_old <- as.vector(
  table(mmus$Identity_ref_fraction[mmus$Age_ID == "old"]))[
    match(celltype_df_dol$celltype, names(table(
      mmus$Identity_ref_fraction[mmus$Age_ID == "old"])))]

celltype_df_dol$mcas_yng <- as.vector(
  table(mcas$Identity_ref_fraction[mcas$Age_ID == "yng"]))[
    match(celltype_df_dol$celltype, names(table(
      mcas$Identity_ref_fraction[mcas$Age_ID == "yng"])))]

celltype_df_dol$mcas_old <- as.vector(
  table(mcas$Identity_ref_fraction[mcas$Age_ID == "old"]))[
    match(celltype_df_dol$celltype, names(table(
      mcas$Identity_ref_fraction[mcas$Age_ID == "old"])))]

celltype_df_dol$mcar_yng <- as.vector(
  table(mcar$Identity_ref_fraction[mcar$Age_ID == "yng"]))[
    match(celltype_df_dol$celltype, names(table(
      mcar$Identity_ref_fraction[mcar$Age_ID == "yng"])))]

celltype_df_dol$mcar_old <- as.vector(
  table(mcar$Identity_ref_fraction[mcar$Age_ID == "old"]))[
    match(celltype_df_dol$celltype, names(table(
      mcar$Identity_ref_fraction[mcar$Age_ID == "old"])))]

celltype_df_dol$mspr_yng <- as.vector(
  table(mspr$Identity_ref_fraction[mspr$Age_ID == "yng"]))[
    match(celltype_df_dol$celltype, names(table(
      mspr$Identity_ref_fraction[mspr$Age_ID == "yng"])))]

celltype_df_dol$mspr_old <- as.vector(
  table(mspr$Identity_ref_fraction[mspr$Age_ID == "old"]))[
    match(celltype_df_dol$celltype, names(table(
      mspr$Identity_ref_fraction[mspr$Age_ID == "old"])))]

celltype_df_dol[is.na(celltype_df_dol)] <- 0

# Dahlin
celltype_df_dahl <- data.frame(
  celltype = names(table(sce_dah$identity_ref))
)

celltype_df_dahl$mmus_yng <- as.vector(
  table(mmus$Identity_ref_fraction[mmus$Age_ID == "yng"]))[
    match(celltype_df_dahl$celltype, names(table(
      mmus$Identity_ref_fraction[mmus$Age_ID == "yng"])))]

celltype_df_dahl$mmus_old <- as.vector(
  table(mmus$Identity_ref_fraction[mmus$Age_ID == "old"]))[
    match(celltype_df_dahl$celltype, names(table(
      mmus$Identity_ref_fraction[mmus$Age_ID == "old"])))]

celltype_df_dahl$mcas_yng <- as.vector(
  table(mcas$Identity_ref_fraction[mcas$Age_ID == "yng"]))[
    match(celltype_df_dahl$celltype, names(table(
      mcas$Identity_ref_fraction[mcas$Age_ID == "yng"])))]

celltype_df_dahl$mcas_old <- as.vector(
  table(mcas$Identity_ref_fraction[mcas$Age_ID == "old"]))[
    match(celltype_df_dahl$celltype, names(table(
      mcas$Identity_ref_fraction[mcas$Age_ID == "old"])))]

celltype_df_dahl$mcar_yng <- as.vector(
  table(mcar$Identity_ref_fraction[mcar$Age_ID == "yng"]))[
    match(celltype_df_dahl$celltype, names(table(
      mcar$Identity_ref_fraction[mcar$Age_ID == "yng"])))]

celltype_df_dahl$mcar_old <- as.vector(
  table(mcar$Identity_ref_fraction[mcar$Age_ID == "old"]))[
    match(celltype_df_dahl$celltype, names(table(
      mcar$Identity_ref_fraction[mcar$Age_ID == "old"])))]

celltype_df_dahl$mspr_yng <- as.vector(
  table(mspr$Identity_ref_fraction[mspr$Age_ID == "yng"]))[
    match(celltype_df_dahl$celltype, names(table(
      mspr$Identity_ref_fraction[mspr$Age_ID == "yng"])))]

celltype_df_dahl$mspr_old <- as.vector(
  table(mspr$Identity_ref_fraction[mspr$Age_ID == "old"]))[
    match(celltype_df_dahl$celltype, names(table(
      mspr$Identity_ref_fraction[mspr$Age_ID == "old"])))]

celltype_df_dahl[is.na(celltype_df_dahl)] <- 0

# RBIND
celltype_df_fraction <- rbind(celltype_df_dol, celltype_df_dahl)

# Baccin 
celltype_df_bac <- data.frame(
  celltype = names(table(sce_bac$identity_ref))
)

celltype_df_bac$mmus_yng <- as.vector(
  table(mmus$Identity_ref_all[mmus$Age_ID == "yng"]))[
    match(celltype_df_bac$celltype, names(table(
      mmus$Identity_ref_all[mmus$Age_ID == "yng"])))]

celltype_df_bac$mmus_old <- as.numeric(
  table(mmus$Identity_ref_all[mmus$Age_ID == "old"]))[
    match(celltype_df_bac$celltype, names(table(
      mmus$Identity_ref_all[mmus$Age_ID == "old"])))]

celltype_df_bac$mcas_yng <- as.numeric(
  table(mcas$Identity_ref_all[mcas$Age_ID == "yng"]))[
    match(celltype_df_bac$celltype, names(table(
      mcas$Identity_ref_all[mcas$Age_ID == "yng"])))]

celltype_df_bac$mcas_old <- as.numeric(
  table(mcas$Identity_ref_all[mcas$Age_ID == "old"]))[
    match(celltype_df_bac$celltype, names(table(
      mcas$Identity_ref_all[mcas$Age_ID == "old"])))]

celltype_df_bac$mcar_yng <- as.numeric(
  table(mcar$Identity_ref_all[mcar$Age_ID == "yng"]))[
    match(celltype_df_bac$celltype, names(table(
      mcar$Identity_ref_all[mcar$Age_ID == "yng"])))]

celltype_df_bac$mcar_old <- as.numeric(
  table(mcar$Identity_ref_all[mcar$Age_ID == "old"]))[
    match(celltype_df_bac$celltype, names(table(
      mcar$Identity_ref_all[mcar$Age_ID == "old"])))]

celltype_df_bac$mspr_yng <- as.numeric(
  table(mspr$Identity_ref_all[mspr$Age_ID == "yng"]))[
    match(celltype_df_bac$celltype, names(table(
      mspr$Identity_ref_all[mspr$Age_ID == "yng"])))]

celltype_df_bac$mspr_old <- as.numeric(
  table(mspr$Identity_ref_all[mspr$Age_ID == "old"]))[
    match(celltype_df_bac$celltype, names(table(
      mspr$Identity_ref_all[mspr$Age_ID == "old"])))]

celltype_df_bac[is.na(celltype_df_bac)] <- 0

library(tidyverse)

celltype_df_bac <- column_to_rownames(celltype_df_bac, var = "celltype")
celltype_df_bac <- rbind(celltype_df_bac, colSums(celltype_df_bac))
rownames(celltype_df_bac)[nrow(celltype_df_bac)] <- "SUM"
celltype_df_bac <- rownames_to_column(celltype_df_bac, var = "celltype")

celltype_df_fraction <- column_to_rownames(celltype_df_fraction, var = "celltype")
celltype_df_fraction <- rbind(celltype_df_fraction, colSums(celltype_df_fraction))
rownames(celltype_df_fraction)[nrow(celltype_df_fraction)] <- "SUM"
celltype_df_fraction <- rownames_to_column(celltype_df_fraction, var = "celltype")

write_csv(celltype_df_bac, "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/celltypenrs_refbaccin.csv")
write_csv(celltype_df_fraction, "/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/celltypenrs_refdolgalev_dahlin.csv")





mspr <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/09_seurat3/sce_mspr_Batch_exp_day-09")
mmus <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/09_seurat3/sce_mmus_Batch_exp_day-09")
mcar <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/09_seurat3/sce_mcar_Batch_exp_day-09")
mcas <- readRDS(file = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/09_seurat3/sce_mcas_Batch_exp_day-09")
mspr <- reduce_dims(mspr, nr_hvgs = 2000)
mmus <- reduce_dims(mmus, nr_hvgs = 2000)
mcar <- reduce_dims(mcar, nr_hvgs = 2000)
mcas <- reduce_dims(mcas, nr_hvgs = 2000)

umap_base(mspr, color_by = "Batch_exp_day")
umap_base(mmus, color_by = "Batch_exp_day")
umap_base(mcar, color_by = "Batch_exp_day")
umap_base(mcas, color_by = "Batch_exp_day")

