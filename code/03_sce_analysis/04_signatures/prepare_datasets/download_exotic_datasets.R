#-------------------------------------------------------------------------------

# download datasets containing bone marrow cells (2 human HSPC, 2 human stromal)

set.seed(37)
options(timeout=1000)

# get input paths
Seurat_hgl_sorted_BM <- snakemake@input[["Seurat_hgl_sorted_BM"]]
Seurat_hgl_whole_BM <- snakemake@input[["Seurat_hgl_whole_BM"]]


#-------------------------------------------------------------------------------
# naked mole rat = NMR = heterocephalus glaber = HGL
# downloaded from: https://figshare.com/collections/Naked_mole-rat_HSPC_Single-cell_RNA-Seq_Seurat_objects_R/5474256
# https://doi.org/10.6084/m9.figshare.c.5474256
# https://doi.org/10.15252%2Fembj.2021109694

utils::download.file(
  url = "https://figshare.com/ndownloader/files/28484826",
  destfile = Seurat_hgl_sorted_BM)

utils::download.file(
  url = "https://figshare.com/ndownloader/files/28484826",
  destfile = Seurat_hgl_sorted_BM)

#-------------------------------------------------------------------------------
# zebrafish = zeb = danio rerio = dnr
# downloaded from: https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-5530/downloads?ref=biostudies
# https://doi.org/10.1038/s41467-017-02305-6

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download?fileType=cluster&accessKey=",
  destfile = zeb_clustering_file)

#zeb_clustering_file <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw/zebrafish/zeb_clusters.tsv"
# E-MTAB-5530.clusters.tsv

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download/zip?fileType=quantification-filtered&accessKey=",
  destfile = zeb_quantification_filtered_files)

#zeb_quantification_filtered_files <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw/zebrafish//zeb_quantification_filtered_files.zip"
# E-MTAB-5530-quantification-filtered-files.zip

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download/zip?fileType=normalised&accessKey=",
  destfile = zeb_normalised_files)

#zeb_normalised_files <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw/zebrafish/zeb_normalised_files.zip"
# E-MTAB-5530-normalised-files.zip

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download/zip?fileType=quantification-raw&accessKey=",
  destfile = zeb_quantification_raw_files)

#zeb_quantification_raw_files <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw/zebrafish/zeb_quantification_raw_files.zip"
# E-MTAB-5530-quantification-raw-files.zip