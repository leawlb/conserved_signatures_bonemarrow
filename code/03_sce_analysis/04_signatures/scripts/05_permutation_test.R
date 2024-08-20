
# perform permutation tests on reclustering scores as chosen


#-------------------------------------------------------------------------------
set.seed(37)

library(Seurat, quietly = TRUE)
library(parallel, quietly = TRUE)
library(dplyr, quietly = TRUE)

set.seed(37)
base::RNGkind("L'Ecuyer-CMRG") # for mclapply, generation of same random numbers

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD OBJECTS

# params for re-clustering and permutation test

# cut_off_counts <- 5
# nr_cores <- 1
# iterations <- 2

cut_off_counts <- snakemake@params[["cut_off_counts"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

#-------------------------------------------------------------------------------

# decide which gene sets to use for the current dataset (based on fraction)
dataset_curr <- snakemake@wildcards[["dataset"]]

# dataset_curr <- "li_all_stromal"

if(dataset_curr %in% c("ts_all_stromal", "li_all_stromal")){
  fraction_curr <- "str"
}else if(dataset_curr %in% c("ts_bone_marrow", 
                             "ts_hscs_progenitors",
                             "nmr_sorted_hspc",
                             "zeb_all_hspc")){
  fraction_curr <- "hsc"
}

#-------------------------------------------------------------------------------

# load the seurat object/dataset to be tested
seu_preprocessed <- base::readRDS(snakemake@input[["seu_preprocessed"]])

#seu_preprocessed <- base::readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/prepared/li_all_stromal")

#-------------------------------------------------------------------------------

# load ensembl IDs of specific gene set to be tested from correct fraction
# this can be signature genes or all BL6 marker genes
ensembl_paths <- snakemake@input[["ensembl_paths"]]
#ensembl_paths <- c("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/02_endf/ensembl_sign_str")
print(ensembl_paths)

# this is the correct one for conserved_signature_genes
ensembl_df <- base::readRDS(ensembl_paths[[
  which(base::grepl(fraction_curr, ensembl_paths))]])
print(head(ensembl_df))

#-------------------------------------------------------------------------------

# load the datafram which contains info on which resolution to use for
# re-clustering
resolution_df_path <- snakemake@params[["resolution_df"]] 
#resolution_df_path <- "/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/reclustering_other_resolution.txt"

resolution_df <- utils::read.csv(file = resolution_df_path, 
                                 header = TRUE, 
                                 sep = ";", 
                                 check.names=FALSE, 
                                 stringsAsFactors=FALSE, 
                                 as.is=TRUE, 
                                 colClasses = "character")

resolution_df <- resolution_df[resolution_df$dataset == dataset_curr,]

resl <- resolution_df$resolution[
  resolution_df$conservation_level == "random_features"]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PREPARE

# extract conserved signature gene set, keep only genes that are also in the 
# Seurat Object 
ensembl_column_use <- seu_preprocessed@misc$ensembl_column_use

conserved_signature_IDs <- base::unique(
  ensembl_df[,which(colnames(ensembl_df) == ensembl_column_use)])
conserved_signature_IDs <- conserved_signature_IDs[
  conserved_signature_IDs %in% rownames(seu_preprocessed)]

# get the number of conserved signature genes, the same number of random genes 
# will be used for the permutation test
nr_recl_genes <- length(conserved_signature_IDs)
print(nr_recl_genes)
print(dataset_curr)

#-------------------------------------------------------------------------------

# get only genes that have a count of at least n = cut_off_counts
print(base::summary(rowSums(seu_preprocessed@assays$RNA$counts)))

gene_pool <- rownames(seu_preprocessed)[
  which(rowSums(seu_preprocessed@assays$RNA$counts) > cut_off_counts)]
print(length(gene_pool))

# subset seurat object to these genes as the other genes won't be used
seu_preprocessed_sub <- SeuratObject::subset(seu_preprocessed,
                                             features = gene_pool,
                                             slot = "count")
print(dim(seu_preprocessed_sub))
print(base::summary(rowSums(seu_preprocessed_sub@assays$RNA$counts)))

#-------------------------------------------------------------------------------

# generate i = iterations random sets of the same length as there are emfs
# always generate the same random numbers
set.seed(37)
base::RNGkind("L'Ecuyer-CMRG")

iteration_df <- base::data.frame(row.names = c(1:nr_recl_genes))

for(i in 1:iterations){
  iteration_df[,i] <- base::sample(1:length(gene_pool), 
                                   nr_recl_genes, 
                                   replace = FALSE)
  colnames(iteration_df)[i] <- i
}
print(iteration_df[1:10,1:3])

#-------------------------------------------------------------------------------
# run standard seurat re-clustering pipeline (own function)

print(nr_cores)
print(cut_off_counts)
print(resl)

res_df_list <- parallel::mclapply(X = as.list(c(1:iterations)),
                                  FUN = random_reclustering_scores,
                                  seu = seu_preprocessed_sub,
                                  assay_use = "RNA",
                                  iteration_df = iteration_df,
                                  resolution = resl,
                                  mc.preschedule = TRUE, 
                                  mc.cores = nr_cores,
                                  mc.silent = TRUE)

# res_df_list <- lapply(X = as.list(c(1:iterations)),
#                                   FUN = random_reclustering_scores,
#                                   seu = seu_preprocessed_sub,
#                                   assay_use = "RNA",
#                                   iteration_df = iteration_df,
#                                   resolution = resolution)

score_df <- dplyr::bind_rows(res_df_list)
print(head(score_df))

#-------------------------------------------------------------------------------

base::saveRDS(score_df, snakemake@output[["perm_score_df"]])

utils::sessionInfo()
