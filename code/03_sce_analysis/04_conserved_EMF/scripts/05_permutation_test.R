
# perform permutation tests on reclustering scores
# score 1: proportion of fields that are 0 from the max possible nr of fields 
#          that can be 0 in a theoreical perfect re-clustering.
# score 2: the mean of the proportion of cell per cell type and cell per cluster
#          for each field in the comparison matrix.
# In a perfect re-clustering, both scores would be 1

#-------------------------------------------------------------------------------
set.seed(37)

library(scCustomize)
library(Seurat)
library(parallel)
library(dplyr)

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------

if(snakemake@wildcards[["reference"]] %in% 
   c("ts_all_stromal", "li_all_stromal")){
  fraction_curr <- "str"
}else if(snakemake@wildcards[["reference"]] %in% 
         c("ts_bone_marrow", "ts_hsc_progenitors")){
  fraction_curr <- "hsc"
}

# load the human reference dataset
seu_preprocessed <- readRDS(snakemake@input[["seu_preprocessed"]])
#seu_preprocessed <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/metadata/scRNAseq/03_sce_analysis/reclustering_bm/pre-processed/ts_all_stromal")
emfs_paths <- snakemake@input[["ensembl_emfs"]]
ensembl_emfs <- readRDS(emfs_paths[[which(!grepl(fraction_curr, emfs_paths))]])
#ensembl_emfs <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_conserved_EMF/03_ensm/ensembl_emfs_str")
nr_emfs <- nrow(ensembl_emfs)
print(snakemake@wildcards[["reference"]] )
print(nr_emfs)

cut_off_counts <- snakemake@params[["cut_off_counts"]]
iterations <- snakemake@params[["iterations"]]
# iterations <- 3
resolution_list <- snakemake@params[["resolution"]] 
resolution <- resolution_list[[snakemake@wildcards[["reference"]]]]
# resolution <- 0.4
nr_cores <- snakemake@params[["nr_cores"]] 
# nr_cores <- 2

#-------------------------------------------------------------------------------
# generate i sets of genes that are expressed in at 
# least 10% of cells of at least 1 cell type
# this is the pre-requisite for marker genes from Seurat::FindMarkers()

# get the percentage of cells per cell type
#perc_expr <- scCustomize::Percent_Expressing(seu_preprocessed, layer = "counts",
#                                             features = Features(seu_preprocessed), 
#                                             split_by = "cell_type")

# get the genes expressed in at least 10% of cells of at least 1 cell type
#gene_pool <- vector()
#for(i in 1:ncol(perc_expr)){
#  gene_pos <- which(perc_expr[,i] > 5)
#  print(length(gene_pos))
#  gene_pool <- c(gene_pool, rownames(perc_expr)[gene_pos])
#}
#gene_pool <- unique(gene_pool)
#print(length(gene_pool))

#-------------------------------------------------------------------------------
# get the nr of genes that have a counts of at least 10
summary(rowSums(seu_preprocessed@assays$RNA$counts))

gene_pool <- rownames(seu_preprocessed)[
  which(rowSums(seu_preprocessed@assays$RNA$counts) > cut_off_counts)]
print(length(gene_pool))

# subset seurat object to these genes as the other genes won't be used
seu_preprocessed_sub <- subset(seu_preprocessed, features = gene_pool,
                               slot = "count")
dim(seu_preprocessed_sub)
summary(rowSums(seu_preprocessed_sub@assays$RNA$counts))

# generate i = iterations random sets of the same length as there are emfs
# always generate the same random numbers
set.seed(37)
RNGkind("L'Ecuyer-CMRG")
iteration_df <- data.frame(row.names = c(1:nr_emfs))
for(i in 1:iterations){
  iteration_df[,i] <- sample(1:length(gene_pool), 
                             nr_emfs, 
                             replace = FALSE)
  colnames(iteration_df)[i] <- i
}
print(iteration_df[1:10,1:3])

#-------------------------------------------------------------------------------

print(nr_cores)
print(cut_off_counts)
print(resolution)

# own function
res_df_list <- parallel::mclapply(X = as.list(c(1:iterations)),
                                  FUN = random_reclustering_scores,
                                  seu = seu_preprocessed_sub,
                                  assay_use = "RNA",
                                  iteration_df = iteration_df,
                                  resolution = resolution,
                                  mc.preschedule = TRUE, 
                                  mc.cores = nr_cores,
                                  mc.silent = TRUE)
score_df <- bind_rows(res_df_list)

saveRDS(score_df, snakemake@output[["perm_score_df"]])
