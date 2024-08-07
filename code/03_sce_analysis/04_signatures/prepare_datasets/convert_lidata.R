#-------------------------------------------------------------------------------

# script to convert data from Li et al. publication to SCE via zellkonverter, 
# then from SCE to Seurat for consistency with other test datasets

#-------------------------------------------------------------------------------
# load libraries, objects

library(zellkonverter, quietly = TRUE)
library(basilisk, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(Seurat, quietly = TRUE)

set.seed(37)

h5ad_input <- snakemake@input[["h5ad_input"]]
h5ad_raw_input <- snakemake@input[["h5ad_raw_input"]]
print(h5ad_input)
print(h5ad_raw_input)

#-------------------------------------------------------------------------------
# conversion from h5ad format/adata object

# convert and load in one step from zellkonverter
sce <- zellkonverter::readH5AD(h5ad_input)
sce_raw <- zellkonverter::readH5AD(h5ad_raw_input)

dim(sce)
dim(sce_raw)

# safekeeping
base::saveRDS(sce, snakemake@output[["sce_output_1"]])
base::saveRDS(sce_raw, snakemake@output[["sce_raw_output_1"]])

#-------------------------------------------------------------------------------
# prepare sce for seurat conversion: assays, names

# remove assays from sce object that are not required
SummarizedExperiment::assays(sce)$Mu <- NULL
SummarizedExperiment::assays(sce)$Ms <- NULL
SummarizedExperiment::assays(sce)$spliced <- NULL
SummarizedExperiment::assays(sce)$unspliced <- NULL

# remove duplicated genes from the sce object
dup_genes <- base::unique(rownames(sce)[which(base::duplicated(rownames(sce)))])
sce <- sce[-which(rownames(sce) %in% dup_genes),]

# add colnames to the object
colnames(sce) <- sce$Unnamed..0

base::table(is.na(base::match(rownames(sce), rownames(sce_raw))))
base::table(is.na(base::match(colnames(sce), colnames(sce_raw))))

# match sce_raw to sce
sce_raw <- sce_raw[match(rownames(sce), rownames(sce_raw)),
                   match(colnames(sce), colnames(sce_raw))]

# transfer assay of raw counts to sce
SummarizedExperiment::assays(sce)$counts <- SummarizedExperiment::assays(sce_raw)$X 
# put in a placeholder for logcounts required for Seurat conversion
SummarizedExperiment::assays(sce)$logcounts <- SummarizedExperiment::assays(sce)$X

# for safekeeping also
base::saveRDS(sce, snakemake@output[["sce_output_2"]])

#-------------------------------------------------------------------------------
# add cell type info from publication

# add info on clusters from publication figures 1 and 2
non_hematopoietic_clusters <- c(39, 0, 5, 29, 6, 38, 28, 16, 23, 37, 8, 3)

sce_nh <- sce[,sce$louvain11_merged.95_renamed %in% non_hematopoietic_clusters]
sce_nh$cell_type_from_publication <- vector(length = ncol(sce_nh))

# detailed info on stromal cell types from publication figure 2
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 3] <- 'multipotent stromal stem cells' # MSSC 
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 5] <- 'highly adipocytic gene-expressing progenitors' # HAGEP
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 16] <- 'balanced prog.'
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 38] <- 'pre-osteoblast'
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 29] <- 'osteochondrogenic progenitors 1' # OC
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 23] <- 'osteochondrogenic progenitors 2' # OC
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 6] <- 'pre-fibroblast 1'
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 37] <- 'pre-fibroblast 2'
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 8] <- 'pre-fibroblast 3'

# info from figure 1
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 28] <- 'endothelial'
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 0] <- 'basal'
sce_nh$cell_type_from_publication[sce_nh$louvain11_merged.95_renamed == 39] <- 'neuronal'

# remove cell types that are not in our dataset
cts_remove <- c('basal', 'neuronal')

sce_nh <- sce_nh[,which(!sce_nh$cell_type_from_publication %in% cts_remove)]
table(sce_nh$cell_type_from_publication)

#-------------------------------------------------------------------------------
# change rownames from gene symbols to IDs as required downstream
rownames(sce_nh) <- rowData(sce_nh)$GeneID

# safekeeping
base::saveRDS(sce_nh, snakemake@output[["sce_output_3"]])

#-------------------------------------------------------------------------------
# conversion to Seurat object for consistency

# convert 
seu_nh <- Seurat::as.Seurat(sce_nh)

# consistency with other test datasets
seu_nh@neighbors <- list()
seu_nh@reductions <- list()

RNA <- Seurat::CreateAssayObject(seu_nh@assays$originalexp$counts)
seu_nh[['RNA']] <- RNA
Seurat::DefaultAssay(seu_nh) <- 'RNA'

seu_nh$cell_type <- seu_nh$cell_type_from_publication

seu_nh$cell_type <- factor(seu_nh$cell_type, 
                           levels = c('multipotent stromal stem cells',
                                      'balanced prog.',
                                      'highly adipocytic gene-expressing progenitors',
                                      'pre-osteoblast',
                                      'osteochondrogenic progenitors 1',
                                      'osteochondrogenic progenitors 2',
                                      'pre-fibroblast 1',
                                      'pre-fibroblast 2',
                                      'pre-fibroblast 3',
                                      'endothelial'))

base::table(seu_nh$cell_type)

base::saveRDS(seu_nh, snakemake@output[["seurat_output"]])

utils::sessionInfo()