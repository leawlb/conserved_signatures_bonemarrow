
#-------------------------------------------------------------------------------

# TABULA MURIS DATASET
# TODO: add to snakefile

# https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687

#-------------------------------------------------------------------------------

url_metadata_all <- "https://figshare.com/ndownloader/files/10842785"
url_annotation_all <- "https://figshare.com/ndownloader/files/13088129"

path_metadata_all <- snakemake@output[["path_metadata_all"]]
path_annotation_all <- snakemake@output[["path_annotation_all"]]

# download metadata
download.file(
  url = url_metadata_all,
  destfile = path_metadata_all)

# download annotation
download.file(
  url = url_annotation_all, 
  destfile = path_annotation_all)

# download FACS folder manually and upload or move into correct folder
# https://figshare.com/ndownloader/files/10700143


#-------------------------------------------------------------------------------

# WEINREB

# GSE140802
# https://pubmed.ncbi.nlm.nih.gov/31974159/

options(timeout=5000)

url_weinreb <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE140802&format=file"
# path_weinreb <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/weinreb/GSE140802_RAW.tar"

path_weinreb <- snakemake@output[["path_weinreb"]]

# this takes long
download.file(url = url_weinreb,
              destfile = path_weinreb)
