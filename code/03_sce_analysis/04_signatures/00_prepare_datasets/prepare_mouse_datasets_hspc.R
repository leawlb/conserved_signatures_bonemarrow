
# prepare mouse datasets

# determine random number generator for sample
# Mersenne-Twister" is default
base::RNGkind("Mersenne-Twister")

set.seed(37)

library(Matrix)
library(Seurat)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# TABULA MURIS DATASETS

# https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687

print("TABULA MURIS DATASETS")

#-------------------------------------------------------------------------------

# load raw counts
#path_counts_marrow <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/FACS/Marrow-counts.csv"
path_counts_marrow <- snakemake@input[["path_counts_marrow"]]
print("counts_marrow")
print(path_counts_marrow)

counts_marrow <- utils::read.csv(
  file = path_counts_marrow, 
  header = TRUE, 
  sep = ",", 
  check.names=FALSE, 
  stringsAsFactors=FALSE, 
  as.is=TRUE)

print(counts_marrow[1:5,1:5])
print(dim(counts_marrow))
print(head(rownames(counts_marrow)))

# put gene names into rownames 
rownames(counts_marrow) <- counts_marrow[,1]
counts_marrow <- counts_marrow[,-1]

# remove genes with underscores to avoid error later (only four genes)
counts_marrow <- counts_marrow[-grep("_", rownames(counts_marrow)),]

print(dim(counts_marrow))

# convert all counts to numerical 
for(c in 1:ncol(counts_marrow)){
  counts_marrow[,c] <- as.numeric(counts_marrow[,c])
}

# convert to sparse matrix
counts_marrow <- base::as.matrix(counts_marrow)
counts_marrow <- as(counts_marrow, "dgCMatrix")
print(counts_marrow[1:10, 1:10])

#-------------------------------------------------------------------------------

# load annotation with cell type annotations
#path_annotation_all <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/annotations_facs.csv"
path_annotation_all <- snakemake@input[["path_annotation_all"]]
print("annotation_all")
print(path_annotation_all)

annotation_all <- utils::read.csv(
  file = path_annotation_all, 
  header = TRUE, 
  sep = ",", 
  check.names=FALSE, 
  stringsAsFactors=FALSE, 
  as.is=TRUE, 
  colClasses = "character")

print(annotation_all[1:5,1:5])
print(dim(annotation_all))
print(colnames(annotation_all))
print(annotation_all[1:5,1:2])
annotation_all <- annotation_all[,-c(1,2)]

print(base::table(!is.na(base::match(colnames(counts_marrow),
                                     annotation_all$cell))))
print(base::table(!is.na(base::match(annotation_all$cell, 
                                     colnames(counts_marrow)))))

# subset annotation to only relevant cell types/tissues
print(base::table(annotation_all$tissue))
annotation_marrow <- annotation_all[annotation_all$tissue == "Marrow",]
print(base::table(annotation_marrow$cell_ontology_class))

# put cell barcodes into rows of annotation
rownames(annotation_marrow) <- annotation_marrow$cell

#-------------------------------------------------------------------------------
# load sample-wise metadata
#path_metadata_all <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/metadata_FACS.csv"
path_metadata_all <- snakemake@input[["path_metadata_all"]]
print("metadata_all")
print(path_metadata_all)

metadata_all <- utils::read.csv(file = path_metadata_all, 
                                  header = TRUE, 
                                  sep = ",", 
                                  check.names=FALSE, 
                                  stringsAsFactors=FALSE, 
                                  as.is=TRUE, 
                                  colClasses = "character")

print(metadata_all[1:4,1:4])
print(dim(metadata_all))

# maybe not required for now but definitely good to have

#-------------------------------------------------------------------------------
# generate seurat object

print("preparing for generating seurat")

print(annotation_marrow[1:4, 1:4])
print(counts_marrow[1:4, 1:4])

print("generating seurat")

seu_tmu <- SeuratObject::CreateSeuratObject(
  meta.data = annotation_marrow,
  counts = counts_marrow, 
  project = "SeuratProject", 
  assay = "RNA")
# no reduced dimensions, no logcounts, only raw counts

print(seu_tmu)
print(seu_tmu@assays)
print(head(seu_tmu@meta.data))

#-------------------------------------------------------------------------------
# prepare seurat object

base::table(seu_tmu@meta.data$cell_ontology_class)

ct_remove <- c(
  "B cell",
  "basophil",
  "late pro-B cell",
  "granulocyte",
  "macrophage",
  "mature natural killer cell",
  "monocyte",
  "naive B cell",
  "pre-natural killer cell",
  "precursor B cell",
  "regulatory T cell",
  "immature natural killer cell",
  "immature NK T cell",
  "immature T cell"
)

seu_tmu <- seu_tmu[,which(!seu_tmu$cell_ontology_class %in% ct_remove)]

base::table(seu_tmu$cell_ontology_class)
base::table(is.na(seu_tmu$cell_ontology_class))

# remove NAs
seu_tmu <- seu_tmu[,which(!seu_tmu$cell_ontology_class == "NA")]

# put cell types into "cell_type" slot
seu_tmu$cell_type <- seu_tmu$cell_ontology_class

seu_tmu$cell_type <- factor(seu_tmu$cell_type, levels = c(
  "Slamf1-positive multipotent progenitor cell",
  "Slamf1-negative multipotent progenitor cell",
  "hematopoietic precursor cell",
  "granulocytopoietic cell",
  "granulocyte monocyte progenitor cell",
  "megakaryocyte-erythroid progenitor cell",
  "common lymphoid progenitor",
  "immature B cell"
))

stopifnot(!is.na(seu_tmu$cell_type))
print(base::table(seu_tmu$cell_type))

#-------------------------------------------------------------------------------

# check rownames
rownames(seu_tmu)[!grepl("Rik", rownames(seu_tmu))][200:230]
# add info on which column of the ensembl data frame to use based on Features 
seu_tmu@misc$ensembl_column_use <- "MMUS_SYMBOL" # mouse gene symbols

#-------------------------------------------------------------------------------
# add info on which assay should be used for reclustering
seu_tmu@misc$data_use <- "raw_counts"

#-------------------------------------------------------------------------------

#base::saveRDS(seu_tmu, "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/prepared/tm_bone_marrow")
base::saveRDS(seu_tmu, snakemake@output[["mus_tm_bonemarrow"]])


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# WEINREB

print("WEINREB")

#-------------------------------------------------------------------------------
# load metadata file with annotations

#path_weinreb_base <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/weinreb/"
path_weinreb_base <- snakemake@input[["path_weinreb_base"]]
print(path_weinreb_base)

print("metadata")
# use in vivo data = hscs that were partly differentiated in recipient mice
metadata <- utils::read.table(
  base::paste0(path_weinreb_base, "/GSM4185643_stateFate_inVivo_metadata.txt.gz"), 
  stringsAsFactors = FALSE,
  header = TRUE,
  sep = "\t")

print(base::table(metadata$Cell.type.annotation))

#-------------------------------------------------------------------------------
# load library names (maybe not required)

print("library_names")
library_names <- utils::read.table(
  base::paste0(path_weinreb_base, "/GSM4185643_stateFate_inVivo_library_names.txt.gz"), 
  stringsAsFactors = FALSE,
  sep = "\t")

print(base::table(library_names[,1] %in% metadata$Library))
print(base::table(metadata$Library %in% library_names[,1]))
print(base::table(library_names[,1]))

#-------------------------------------------------------------------------------
# load normalized counts

print("normed_counts_matrix")
normed_counts_matrix <- Seurat::ReadMtx(
  mtx = base::paste0(path_weinreb_base, "/GSM4185643_stateFate_inVivo_normed_counts.mtx.gz"), 
  cells = base::paste0(path_weinreb_base, "/GSM4185643_stateFate_inVivo_cell_barcodes.txt.gz"),
  features = paste0(path_weinreb_base, "/GSM4185643_stateFate_inVivo_gene_names.txt.gz"),
  feature.column = 1,
  mtx.transpose = TRUE)

print(normed_counts_matrix[1:10, 1:10])
# normalized counts

#-------------------------------------------------------------------------------
# there are duplicated barcodes
# make them unique in metadata by pasting the Library to each barcode

print(base::table(base::duplicated(metadata$Cell.barcode)))

metadata$unique_barcode <- base::paste0(metadata$Cell.barcode,
                                        "_",
                                        metadata$Library)

print(base::table(base::duplicated(metadata$unique_barcode)))
# 34 duplicated remain

#-------------------------------------------------------------------------------

# make sure that the order of metadata and normed_counts_matrix is the same
base::identical(metadata$Cell.barcode, 
                colnames(normed_counts_matrix))

# get the position of remaining duplicates
dup_pos <- which(base::duplicated((metadata$unique_barcode)))

# remove the positions from both objects
metadata <- metadata[-dup_pos,]
normed_counts_matrix <- normed_counts_matrix[,-dup_pos]

nrow(metadata)
ncol(normed_counts_matrix)

# add unique barcodes as colnames for normed_counts_matrix and as rownames
# for metadata
rownames(metadata) <- metadata$unique_barcode
colnames(normed_counts_matrix) <- metadata$unique_barcode

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# create Seurat object

print("generating seurat object")
seu_wrb <- SeuratObject::CreateSeuratObject(
  meta.data = metadata,
  counts = normed_counts_matrix, 
  project = "SeuratProject", 
  assay = "RNA")
# no reduced dimensions, counts are NORMALIZED!

print(seu_wrb)
print(head(seu_wrb@assays))
print(head(seu_wrb@meta.data))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# prepare

# not all cells are progenitors, some were differentiated in a recipient mouse
print(base::table(seu_wrb$Cell.type.annotation))

# remove cells not in our dataset
ct_remove <- c(
  "T",
  "NK",
  "DC"
)

seu_wrb <- seu_wrb[,which(!seu_wrb$Cell.type.annotation %in% ct_remove)]
print(base::table(seu_wrb$Cell.type.annotation))

# put cell types into "cell_type" slot
seu_wrb$cell_type <- seu_wrb$Cell.type.annotation

seu_wrb$cell_type <- factor(seu_wrb$cell_type, levels = c(
  "Undifferentiated",
  "Prog",
  "Neu",
  "Mono",
  "B",
  "Baso",
  "Ery"
))

print(base::table(seu_wrb$cell_type))
stopifnot(!is.na(seu_wrb$cell_type))

#-------------------------------------------------------------------------------

# check rownames
rownames(seu_wrb)[!grepl("Rik", rownames(seu_wrb))][200:230]

# add info on which column of the ensembl data frame to use based on Features 
seu_wrb@misc$ensembl_column_use <- "MMUS_SYMBOL" # mouse gene symbols


#-------------------------------------------------------------------------------
# add info on which assay should be used for reclustering
seu_wrb@misc$data_use <- "logcounts" # only has normalised counts, so start from there

#-------------------------------------------------------------------------------

# subset to 30,000 random cells
set.seed(777)
subset_pos <- base::sample(c(1:ncol(seu_wrb)), 25000, replace = FALSE)
print(head(subset_pos))
set.seed(37)

seu_wrb <- seu_wrb[,subset_pos]
print(dim(seu_wrb))

#-------------------------------------------------------------------------------

#base::saveRDS(seu_wrb, "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/prepared/weinreb_hspc")
base::saveRDS(seu_wrb, snakemake@output[["mus_weinreb_hspc"]])

utils::sessionInfo()