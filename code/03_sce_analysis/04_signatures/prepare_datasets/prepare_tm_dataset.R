
# tabula muris scRNAseq data

library(Matrix)
library(Seurat)

# https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687
# TODO: add to snakefile
#-------------------------------------------------------------------------------

# raw counts
marrow_counts_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/FACS/Marrow-counts.csv"

marrow_counts <- utils::read.csv(file = marrow_counts_path, 
                                    header = TRUE, 
                                    sep = ",", 
                                    check.names=FALSE, 
                                    stringsAsFactors=FALSE, 
                                    as.is=TRUE)

print(marrow_counts[1:10,1:10])
print(dim(marrow_counts))
print(head(rownames(marrow_counts)))

rownames(marrow_counts) <- marrow_counts[,1]
marrow_counts <- marrow_counts[,-1]

# remove genes with underscores to avoid error (only four genes)
marrow_counts <- marrow_counts[-grep("_", rownames(marrow_counts)),]

print(dim(marrow_counts))

# convert to numerical 
for(c in 1:ncol(marrow_counts)){
  marrow_counts[,c] <- as.numeric(marrow_counts[,c])
}

marrow_counts <- as.matrix(marrow_counts)
marrow_counts <- as(marrow_counts, "dgCMatrix")


#-------------------------------------------------------------------------------

# annotation
annotation_all_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/annotations_facs.csv"

annotation_all <- utils::read.csv(file = annotation_all_path, 
                                 header = TRUE, 
                                 sep = ",", 
                                 check.names=FALSE, 
                                 stringsAsFactors=FALSE, 
                                 as.is=TRUE, 
                                 colClasses = "character")

print(annotation_all[1:10,1:10])
print(dim(annotation_all))

print(table(!is.na(match(colnames(marrow_counts), annotation_all$cell))))
print(table(!is.na(match(annotation_all$cell, colnames(marrow_counts)))))

# dataset has 5037 cells

table(annotation_all$tissue)
annotation_marrow <- annotation_all[annotation_all$tissue == "Marrow",]
table(annotation_marrow$cell_ontology_class)

rownames(annotation_marrow) <- annotation_marrow$cell

#-------------------------------------------------------------------------------
# metadata
metadata_all_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/raw_mus/tabula_muris/metadata_FACS.csv"

metadata_all <- utils::read.csv(file = metadata_all_path, 
                                  header = TRUE, 
                                  sep = ",", 
                                  check.names=FALSE, 
                                  stringsAsFactors=FALSE, 
                                  as.is=TRUE, 
                                  colClasses = "character")

print(metadata_all[1:4,1:4])
print(dim(metadata_all))

# maybe not required for now

#-------------------------------------------------------------------------------
# generate seurat object

seu_tmu <- SeuratObject::CreateSeuratObject(
  meta.data = annotation_marrow,
  counts = marrow_counts, 
  project = "SeuratProject", 
  assay = "RNA")
# no reduced dimensions, no logcounts


head(seu_tmu@assays$RNA$counts)
head(seu_tmu@meta.data)

#-------------------------------------------------------------------------------
# prepare seurat object

base::table(seu_tmu$cell_ontology_class)

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

base::table(seu_tmu$cell_type)

# check rownames
rownames(seu_tmu)[!grepl("Rik", rownames(seu_tmu))][200:230]

# add info on which column of the ensembl data frame to use based on Features 
seu_tmu@misc$ensembl_column_use <- "MMUS_SYMBOL" # mouse gene symbols


saveRDS(seu_tmu, "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/prepared/tm_bone_marrow")
