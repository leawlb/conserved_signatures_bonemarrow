#-------------------------------------------------------------------------------

# download datasets containing bone marrow cells (2 human HSPC, 2 human stromal)

set.seed(37)
options(timeout=1000)

#-------------------------------------------------------------------------------
# get paths from snakemake

tabula_sapiens_bone_marrow <- snakemake@output[["tabula_sapiens_bone_marrow"]]
tabula_sapiens_hsc_progenitors <- snakemake@output[["tabula_sapiens_hsc_progenitors"]]
tabula_sapiens_stromal <- snakemake@output[["tabula_sapiens_stromal"]]

# li is not part of snakemake output, but files are still saved from this script
li_bone_marrow <- snakemake@params[["li_bone_marrow"]]

#-------------------------------------------------------------------------------

# this is the "Bone_Marrow" dataset from tabula sapiens 
# (https://cellxgene.cziscience.com/datasets)
# DOI: 10.1126/science.abl4896
utils::download.file(
  url = "https://datasets.cellxgene.cziscience.com/5e736dcd-01d8-4639-805a-31fea1528be0.rds",
  destfile = tabula_sapiens_bone_marrow)

# this is the "HSC/progenitors" dataset from tabula sapiens
# (https://cellxgene.cziscience.com/datasets)
# DOI: 10.1126/science.abl4896
utils::download.file(
  url = "https://datasets.cellxgene.cziscience.com/367ea9f3-5bc5-4519-bbb4-8c1283e83d71.rds",
  destfile = tabula_sapiens_hsc_progenitors)

# this is the "Stromal cells (all non-immune cells)" dataset from tabula sapiens 
# (https://cellxgene.cziscience.com/datasets)
# DOI: 10.1126/science.abl4896
utils::download.file(
  url = "https://datasets.cellxgene.cziscience.com/085c23b0-6e29-4668-8a19-1673b93fe5b0.rds",
  destfile = tabula_sapiens_stromal)

#-------------------------------------------------------------------------------

# this is a human bone marrow dataset published by Li H et al (2023)
# DOI: https://doi.org/10.7554/elife.81656
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190965
# since these are so many files, they are NOT part of snakemake output
utils::download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Fbarcodes%2Etsv%2Egz",
  destfile = base::paste0(li_bone_marrow, "/download/GSE190965_barcodes.tsv.gz"))
utils::download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Ffeatures%2Etsv%2Egz",
  destfile = base::paste0(li_bone_marrow, "/download/GSE190965_features.tsv.gz"))
utils::download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Fspliced%2Emtx%2Egz",
  destfile = base::paste0(li_bone_marrow, "/download/GSE190965_spliced.mtx.gz"))
utils::download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE190965&format=file&file=GSE190965%5Funspliced%2Emtx%2Egz",
  destfile = base::paste0(li_bone_marrow, "/download/GSE190965_unspliced.mtx.gz"))

utils::download.file(
  url = "https://raw.githubusercontent.com/Hongzhe2022/MSC_BM_scripts/master/data/Hongzhe_clustering_and_2dUMAP_louvain11_merged_95_renamed.csv",
  destfile = base::paste0(li_bone_marrow, "/data/Hongzhe_clustering_and_2dUMAP_louvain11_merged_95_renamed.csv"))

# github repository:
# https://github.com/Hongzhe2022/MSC_BM_scripts/blob/master/scripts/Load_GEO_data_analyze_it.ipynb

# also:
# - manually downloaded Load_GEO_data_analyze_it.ipynb 2024-04-05 https://github.com/Hongzhe2022/MSC_BM_scripts/commit/0e06e526aac2ebb8ca276d3d165f1f7573e7cdbb
# - made a data/ directory in line with github repository structure containing a copy of Hongzhe_clustering_and_2dUMAP_louvain11_merged_95_renamed.csv
# - exported the micromamba environment required to run the jupyter notebook
