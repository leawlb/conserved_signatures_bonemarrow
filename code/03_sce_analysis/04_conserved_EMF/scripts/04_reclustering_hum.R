#-------------------------------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(parallel)

set.seed(37)
RNGkind("L'Ecuyer-CMRG") # for mclapply, generation of same random numbers

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------

cut_off_counts <- snakemake@params[["cut_off_counts"]]
nr_cores <- snakemake@params[["nr_cores"]]

if(snakemake@wildcards[["reference"]] %in% 
   c("ts_all_stromal", "li_all_stromal")){
  fraction_curr <- "str"
}else if(snakemake@wildcards[["reference"]] %in% 
         c("ts_bone_marrow", "ts_hsc_progenitors")){
  fraction_curr <- "hsc"
}
print(snakemake@wildcards[["reference"]])
print(fraction_curr)

seu <- readRDS(snakemake@input[["seu_input"]])
#seu <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/pre-processed/ts_all_stromal")
emfs_paths <- snakemake@input[["ensembl_emfs"]]
mark_paths <- snakemake@input[["ensembl_mark"]]

ensembl_emfs_cor <- readRDS(emfs_paths[[grep(fraction_curr, emfs_paths)]])
#ensembl_emfs_cor <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_conserved_EMF/03_ensm/ensembl_emfs_hsc")
ensembl_emfs_rev <- readRDS(emfs_paths[[which(!grepl(fraction_curr, emfs_paths))]])
#ensembl_emfs_rev <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_conserved_EMF/03_ensm/ensembl_emfs_str")
ensembl_mark_cor <- readRDS(mark_paths[[grep(fraction_curr, mark_paths)]])
#ensembl_mark_cor <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_conserved_EMF/03_ensm/ensembl_mark_hsc")

#-------------------------------------------------------------------------------

# subset by genes

# correct set of emfs
# by using slot = "count" raw counts will be kept in RNA slot 
seu_emf_cor <- subset(seu, features = ensembl_emfs_cor$ensembl_gene_id,
                      slot = "count")

# reverse set of emfs (different fraction)
seu_emf_rev <- subset(seu, features = ensembl_emfs_rev$ensembl_gene_id,
                      slot = "count")

# correct set of marker genes 
seu_mark_cor <- subset(seu, features = ensembl_mark_cor$ensembl_gene_id,
                       slot = "count")

# random set of genes with non-zero expression, same nr as emfs
non_zero_features <- rownames(seu@assays$RNA)[rowSums(seu@assays$RNA) > cut_off_counts]
random_features <- non_zero_features[
  sample(1:length(non_zero_features), 
         length(ensembl_emfs_cor$ensembl_gene_id), replace = F)]
seu_random <- subset(seu, features = random_features,
                     slot = "count")

dim(seu_emf_cor)
dim(seu_emf_rev)
dim(seu_mark_cor)
dim(seu_random)

#-------------------------------------------------------------------------------

# go through *standard* Seurat pipeline 
# made sure that I start from raw counts again 

# to choose optimal cluster resolution later
resolution_vec <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
                    0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.2, 1.3)

print(nr_cores)
print(resolution_vec)
print(cut_off_counts)

#nr_cores <- 1
# run the standard pipeline for each resolution and each subset dataset                           
seu_emf_cor_list <- mclapply(X = as.list(resolution_vec),
                             FUN = standard_seu_pipeline,
                             seu = seu_emf_cor,
                             assay_use = "RNA",
                             calc_umap = TRUE,
                             features = ensembl_emfs_cor$ensembl_gene_id,
                             mc.preschedule = TRUE, 
                             mc.cores = nr_cores,
                             mc.silent = TRUE)
names(seu_emf_cor_list) <- as.character(resolution_vec)

seu_emf_rev_list <- mclapply(X = as.list(resolution_vec),
                           FUN = standard_seu_pipeline,
                           seu = seu_emf_rev,
                           assay_use = "RNA",
                           calc_umap = TRUE,
                           features = ensembl_emfs_rev$ensembl_gene_id,
                           mc.preschedule = TRUE, 
                           mc.cores = nr_cores,
                           mc.silent = TRUE)
names(seu_emf_rev_list) <- as.character(resolution_vec)

seu_mark_cor_list <- mclapply(X = as.list(resolution_vec),
                           FUN = standard_seu_pipeline,
                           seu = seu_mark_cor,
                           assay_use = "RNA",
                           calc_umap = TRUE,
                           features = ensembl_mark_cor$ensembl_gene_id,
                           mc.preschedule = TRUE, 
                           mc.cores = nr_cores,
                           mc.silent = TRUE)
names(seu_mark_cor_list) <- as.character(resolution_vec)

seu_random_list <- mclapply(X = as.list(resolution_vec),
                           FUN = standard_seu_pipeline,
                           seu = seu_random,
                           assay_use = "RNA",
                           calc_umap = TRUE,
                           features = random_features,
                           mc.preschedule = TRUE, 
                           mc.cores = nr_cores,
                           mc.silent = TRUE)
names(seu_random_list) <- as.character(resolution_vec)

export_list <- list("seu_emf_cor" = seu_emf_cor_list,
                    "seu_emf_rev" = seu_emf_rev_list,
                    "seu_mark_cor" = seu_mark_cor_list,
                    "seu_random" = seu_random_list)

saveRDS(export_list, snakemake@output[["seu_output"]])
