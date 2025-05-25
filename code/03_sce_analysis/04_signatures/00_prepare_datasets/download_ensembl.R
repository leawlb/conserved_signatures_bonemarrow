
# downloading ensembl gene ID tables to convert IDs between species

#-------------------------------------------------------------------------------

library(biomaRt, quietly = TRUE)
set.seed(37)

sessionInfo()
packageVersion('biomaRt')
options(timeout=10000)

#-------------------------------------------------------------------------------

# get ensembl mouse and human ID conversion table (download)

# Ensembl release 113 - October 2024
mart_m <- biomaRt::useEnsembl(
  "ensembl",
  dataset="mmusculus_gene_ensembl",
  host="http://oct2024.archive.ensembl.org")

mart_h <- biomaRt::useEnsembl(
  "ensembl",
  dataset="hsapiens_gene_ensembl",
  host="http://oct2024.archive.ensembl.org")

mart_z <- biomaRt::useEnsembl(
  "ensembl",
  dataset="drerio_gene_ensembl",
  host="http://oct2024.archive.ensembl.org")

# keep NMR just for the sake of retaining the same ENSEMBL tables

mart_n <- biomaRt::useEnsembl(
  "ensembl",
  dataset="hgfemale_gene_ensembl",
  host="http://oct2024.archive.ensembl.org")


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

#-------------------------------------------------------------------------------

base::saveRDS(ensembl_mus, snakemake@output[["ensembl_mus"]])
base::saveRDS(ensembl_hum, snakemake@output[["ensembl_hum"]])
base::saveRDS(ensembl_zeb, snakemake@output[["ensembl_zeb"]])
base::saveRDS(ensembl_nmr, snakemake@output[["ensembl_nmr"]])

utils::sessionInfo()
