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
# zebrafish