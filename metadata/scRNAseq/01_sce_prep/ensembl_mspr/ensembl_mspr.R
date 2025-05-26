# This script is required for the download of the ensembl spretus genome info.
# Because spretus SCE objects mapped with the spretus-specific genome don't
# contain rownames.

# Run this from snkmk_isbm environment.
# Last run on 2025-04-16.

#-------------------------------------------------------------------------------

library(biomaRt, quietly = TRUE)

base_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008"

#-------------------------------------------------------------------------------

# download from archived mspr genome Dec 2023 (2023-12-13)
# with help from Perrine Lacour
ensembl_mspr <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                 host="https://jul2023.archive.ensembl.org",
                                 dataset="mspretus_gene_ensembl")

ensembl_list_mspr <- biomaRt::getBM(
  attributes=c("ensembl_gene_id",
               "mmusculus_homolog_ensembl_gene",
               "mmusculus_homolog_associated_gene_name"),
  mart=ensembl_mspr)

saveRDS(ensembl_list_mspr, file = paste0(base_path, "/data/metadata/scRNAseq/01_sce_prep/ensembl_mspr/ensembl_list_mspr"))
