#-------------------------------------------------------------------------------
# Obtain and prepare the reference datasets for cell type annotation 
# more infos in readme.txt

# run this from snkmk_isbm environment (snkmk_isbm.yaml)

#-------------------------------------------------------------------------------

library(Seurat, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(biomaRt, quietly = TRUE)

base_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008"

#-------------------------------------------------------------------------------
# Dahlin et al.
# downloaded 2022-10-28 ABC portal http://abc.sklehabc.com/unicellular/search
metainfo <- read.delim(paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_raw/metaInfo.txt"))
ref_dahlin_seurat <- readRDS(file = paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_raw/notlabel.RDS"))

ref_dahlin_sce <- as.SingleCellExperiment(ref_dahlin_seurat)
ref_dahlin_sce$identity_ref <- metainfo[match(metainfo[,1],
                                              colnames(ref_dahlin_sce)), 2]

unique(ref_dahlin_sce$identity_ref)

# add gene names to ref_dahlin_sce from ensembl list
# use archived version form oct2022
ensembl_mouse <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                         host="https://oct2022.archive.ensembl.org",
                         dataset="mmusculus_gene_ensembl")
ensembl_list_mbl6 <- getBM(attributes=c("ensembl_gene_id",
                                        "external_gene_name"),
                           mart=ensembl_mouse)

rowData(ref_dahlin_sce)$Symbol <- ensembl_list_mbl6$external_gene_name[
  match(rownames(ref_dahlin_sce), ensembl_list_mbl6$ensembl_gene_id)
]

rowData(ref_dahlin_sce)$Symbol[rowData(ref_dahlin_sce)$Symbol == " "] <- NA
rownames(ref_dahlin_sce)[!is.na(rowData(ref_dahlin_sce)$Symbol)] <- rowData(ref_dahlin_sce)$Symbol[!is.na(rowData(ref_dahlin_sce)$Symbol)]

saveRDS(ref_dahlin_sce, file = paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_processed/ref_dahlin_sce"))
sce_test <- readRDS(file =  paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_processed/ref_dahlin_sce"))
unique(sce_test$identity_ref)

#-------------------------------------------------------------------------------
# Dolgalev and Tikhonova
# downloaded 2022-10-28 https://osf.io/ne9vj
ref_dolgalev_seurat <- readRDS(file = paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_raw/bone-marrow-seurat.rds"))
ref_dolgalev_seurat[['RNA']] <- NULL
ref_dolgalev_sce <- as.SingleCellExperiment(ref_dolgalev_seurat)

ref_dolgalev_sce$identity_ref <- ref_dolgalev_sce$label
unique(ref_dolgalev_sce$identity_ref)

saveRDS(ref_dolgalev_sce, file = paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_processed/ref_dolgalev_sce"))
sce_test <- readRDS(file =  paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_processed/ref_dolgalev_sce"))
unique(sce_test$identity_ref)

#-------------------------------------------------------------------------------
# Baccin, Al-Sabah, Velten et al.
# downloaded 2022-10-26 https://nicheview.shiny.embl.de
load(paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_raw/NicheData10x.rda"))
ref_baccin_seurat <- NicheData10x
ref_baccin_sce <- as.SingleCellExperiment(ref_baccin_seurat)

unique(ref_baccin_sce$ident)

ref_baccin_sce$identity_ref <- ref_baccin_sce$ident

saveRDS(ref_baccin_sce, file = paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_processed/ref_baccin_sce"))
sce_test <- readRDS(file = paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/references_processed/ref_baccin_sce"))
unique(sce_test$identity_ref)

