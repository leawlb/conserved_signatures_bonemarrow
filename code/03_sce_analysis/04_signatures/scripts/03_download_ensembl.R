
# downloading ensembl gene ID tables to convert IDs between species

#-------------------------------------------------------------------------------

library(biomaRt, quietly = TRUE)
set.seed(37)
options(timeout=300)

#base_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008"

#-------------------------------------------------------------------------------

# get ensembl mouse and human ID conversion table (download)

mart_m <- biomaRt::useMart(
  "ensembl",
  dataset="mmusculus_gene_ensembl",
  host="https://www.ensembl.org")

mart_h <- biomaRt::useMart(
  "ensembl",
  dataset="hsapiens_gene_ensembl",
  host="https://www.ensembl.org")

mart_z <- biomaRt::useMart(
  "ensembl",
  dataset="drerio_gene_ensembl",
  host="https://www.ensembl.org")

mart_n <- biomaRt::useMart(
  "ensembl",
  dataset="hgfemale_gene_ensembl",
  host="https://www.ensembl.org")

# get ensembl dfs with required attributes

ensembl_mus <- biomaRt::getBM(
  mart = mart_m, 
  attributes = c("ensembl_gene_id", 
                 "external_gene_name"))

ensembl_hum <- biomaRt::getBM(
  mart = mart_h,
  attributes = c("ensembl_gene_id",
                 "external_gene_name", 
                 "mmusculus_homolog_ensembl_gene"))

ensembl_zeb <- biomaRt::getBM(
  mart = mart_z, attributes = c("ensembl_gene_id",
                                "external_gene_name", 
                                "mmusculus_homolog_ensembl_gene"))

ensembl_nmr <- biomaRt::getBM(
  mart = mart_n,
  attributes = c("ensembl_gene_id",
                 "external_gene_name", 
                 "mmusculus_homolog_ensembl_gene"))

print("done")

# base::saveRDS(ensembl_hum, base::paste0(base_path, "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/ensembl_hum"))
# base::saveRDS(ensembl_mus, base::paste0(base_path, "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/ensembl_mus"))
# base::saveRDS(ensembl_zeb, base::paste0(base_path, "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/ensembl_zeb"))
# base::saveRDS(ensembl_nmr, base::paste0(base_path, "/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/ensembl_nmr"))

base::saveRDS(ensembl_mus, snakemake@output[["ensembl_mus"]])
base::saveRDS(ensembl_hum, snakemake@output[["ensembl_hum"]])
base::saveRDS(ensembl_zeb, snakemake@output[["ensembl_zeb"]])
base::saveRDS(ensembl_nmr, snakemake@output[["ensembl_nmr"]])
