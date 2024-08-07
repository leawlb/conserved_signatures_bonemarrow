#-------------------------------------------------------------------------------

# download datasets containing hematopoietic cells from naked mole rat (NMR)
# and zebrafish (ZEB)

set.seed(37)
options(timeout=1000)

# get output paths
Seurat_hgl_sorted_BM <- snakemake@output[["Seurat_hgl_sorted_BM"]]
Seurat_hgl_whole_BM <- snakemake@output[["Seurat_hgl_whole_BM"]]

zeb_clustering_file <- snakemake@output[["zeb_clustering_file"]]
expression_tpm <- snakemake@output[["expression_tpm"]]
aggregated_filtered_normalised_counts <- snakemake@output[["aggregated_filtered_normalised_counts"]]
aggregated_filtered_counts <- snakemake@output[["aggregated_filtered_counts"]]

print(Seurat_hgl_sorted_BM)
print(Seurat_hgl_whole_BM)

print(zeb_clustering_file)
print(expression_tpm)
print(aggregated_filtered_normalised_counts)
print(aggregated_filtered_counts)

#-------------------------------------------------------------------------------
# naked mole rat = NMR = heterocephalus glaber = HGL
# downloaded from: https://figshare.com/collections/Naked_mole-rat_HSPC_Single-cell_RNA-Seq_Seurat_objects_R/5474256
# https://doi.org/10.6084/m9.figshare.c.5474256
# https://doi.org/10.15252%2Fembj.2021109694

utils::download.file(
  url = "https://figshare.com/ndownloader/files/28484826",
  destfile = Seurat_hgl_sorted_BM)

utils::download.file(
  url = "https://figshare.com/ndownloader/files/28484889",
  destfile = Seurat_hgl_whole_BM)

#-------------------------------------------------------------------------------
# zebrafish = zeb = danio rerio = dnr
# downloaded from: https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-5530/downloads?ref=biostudies
# https://doi.org/10.1038/s41467-017-02305-6

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download?fileType=cluster&accessKey=",
  destfile = zeb_clustering_file)
# E-MTAB-5530.clusters.tsv

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download/zip?fileType=quantification-filtered&accessKey=",
  destfile = expression_tpm)
# E-MTAB-5530-quantification-filtered-files.zip
# also called "Filtered TPMs Files" on the website
# TPM = transcript per million, apparently only used for QC 

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download/zip?fileType=normalised&accessKey=",
  destfile = aggregated_filtered_normalised_counts)
# E-MTAB-5530-normalised-files.zip
# also called Normalised Counts File" on the website

utils::download.file(
  url = "https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-5530/download/zip?fileType=quantification-raw&accessKey=",
  destfile = aggregated_filtered_counts)
# E-MTAB-5530-quantification-raw-files.zip
# also called "Raw Counts Files" on the website

utils::sessionInfo()