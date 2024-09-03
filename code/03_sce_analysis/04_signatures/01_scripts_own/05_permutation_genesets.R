
# permutation to compare genes sets (conserved markers and all BL6 markers)
# vs. conserved signature + random genes (same nr) to see if they contain
# important genes in addition to containing conserved signature genes.

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
set.seed(37)
RNGkind("L'Ecuyer-CMRG") 

#-------------------------------------------------------------------------------

library(dplyr, quietly = TRUE)
library(BiocGenerics, quietly = TRUE)
library(S4Vectors, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)

# for using mclapply
library(parallel, quietly = TRUE)

# for calculating scores
library(mclust, quietly = TRUE)
library(mcclust, quietly = TRUE)
library(bluster, quietly = TRUE)
library(dendextend, quietly = TRUE)

source(snakemake@params[["functions_reclustering"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD OBJECTS

# params for re-clustering and permutation test
cut_off_counts <- snakemake@params[["cut_off_counts"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

#-------------------------------------------------------------------------------

# load the seurat object/dataset to be tested
sce_input <- base::readRDS(snakemake@input[["sce_input"]])

fraction_curr <- snakemake@wildcards[["fraction"]]

# cell types to be excluded
cts_exclude <- snakemake@params[["cts_exclude"]]
print(cts_exclude)

# remove cell types to be excluded
sce_input <- sce_input[,which(!sce_input$celltypes %in% cts_exclude)]

#-------------------------------------------------------------------------------
# params for re-clustering

k_graph_list <- snakemake@params[["k_graph_list"]]
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]
resolution_louvain <- resolution_louvain_list[[fraction_curr]]

print(k_graph)
print(resolution_louvain)

#-------------------------------------------------------------------------------
# load geneset list
geneset_list <- base::readRDS(snakemake@input[["geneset_list"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

################   Prepare signature IDs and sce object    #####################

#-------------------------------------------------------------------------------

print("preparing signature IDs")

# get all signature genes, removing subclustering genes in ct specific way
conserved_signature_list <- lapply(geneset_list, function(geneset){
  conserved_signature <- geneset$conserved_signature
  conserved_signature <- conserved_signature[
    which(!conserved_signature %in% geneset$subclustering_genes)]
  return(conserved_signature)
})
sign_IDs <- base::unique(unlist(conserved_signature_list))

nr_sign <- length(sign_IDs) 
print(nr_sign)

#-------------------------------------------------------------------------------
# get position of conserved signature genes 
# POSITIONS FROM ORIGINAL SCE

SIGN_POSITIONS <- which(rownames(sce_input) %in% sign_IDs)
print(length(SIGN_POSITIONS))

#-------------------------------------------------------------------------------
# generate a fitting data frame to add signature gene positions to the 
# random positions downstream

add_df_sign <- base::data.frame(row.names = c(1:nr_sign))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SCE OBJECT
  add_df_sign[,i] <- SIGN_POSITIONS
  colnames(add_df_sign)[i] <- i
}
print("sign_df")
print(add_df_sign[1:10,1:2])
print(dim(add_df_sign))

#-------------------------------------------------------------------------------
# extract subclustering genes

subclustering_gene_list <- lapply(geneset_list, function(geneset){
  subclustering_genes <- geneset$genes_subclustering
  return(subclustering_genes)
})
subclustering_genes <- base::unique(unlist(subclustering_gene_list))

print(length(subclustering_genes))

SUBCLUSTERING_POSITIONS <- which(rownames(sce_input) %in% subclustering_genes)
print(length(SUBCLUSTERING_POSITIONS))

#  check how many subclustering genes became signature genes
# positions in ORIGINAL SCE
print(base::table(SUBCLUSTERING_POSITIONS %in% SIGN_POSITIONS))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

############   Conserved markers vs. signature + random     ####################

#-------------------------------------------------------------------------------

print("Conserved markers vs. signature + random")

# extract conserved marker IDs, removing subclustering genes in ct specific way
conserved_markers_list <- lapply(geneset_list, function(geneset){
  conserved_markers <- geneset$conserved_markers
  conserved_markers <- conserved_markers[
    which(!conserved_markers %in% geneset$subclustering_genes)]
  return(conserved_markers)
})
mark_IDs <- base::unique(unlist(conserved_markers_list))

# get the number of random genes that need to be added to 
# conserved_signature_IDs for each gene set to be permuted
nr_random_mark <- length(mark_IDs) - nr_sign

print(length(mark_IDs))
print(nr_random_mark)

#-------------------------------------------------------------------------------

# subset a SCE object without signature IDs, which contains the pool of 
# genes from which random genes can be drawn:
# - no signature genes
# - no subclustering genes
# - genes expressed n >= cut_off_counts times 

print(length(SIGN_POSITIONS))
print(length(SUBCLUSTERING_POSITIONS))
REMOVE_POS <- base::unique(c(SIGN_POSITIONS, SUBCLUSTERING_POSITIONS))
print(length(REMOVE_POS))
non_sign_sce <- rownames(sce_input)[-REMOVE_POS]

# subset genes
sce_pool <- sce_input[which(rownames(sce_input) %in% non_sign_sce),]
print(nrow(sce_input))
print(nrow(sce_pool))

# get only genes that have a count of at least n = cut_off_counts in RAW counts
gene_pool <- rownames(sce_pool)[
  which(rowSums(assays(sce_pool)$counts) >= cut_off_counts)]
print(length(gene_pool))

# subset seurat object to these genes 
sce_pool <- sce_pool[which(rownames(sce_pool) %in% gene_pool),]

print(dim(sce_pool))
  
#-------------------------------------------------------------------------------

# the positions of genes in SCE that are allowed as pool for random genes
# POSITIONS FROM ORIGINAL SCE
POOL_POSITIONS <- which(rownames(sce_input) %in% rownames(sce_pool))
print(length(POOL_POSITIONS))

#-------------------------------------------------------------------------------

# generate i = iterations random sets of genes at the required number
# from the allowed pool of random genes (original positions in SCE)
# always generate the same random numbers

iteration_df_mark <- base::data.frame(row.names = c(1:nr_random_mark))

set.seed(37)
base::RNGkind("L'Ecuyer-CMRG")
for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SCE subsetted to a smaller list 
  # basically, subset the vector of pool positions in random positions 
  iteration_df_mark[,i] <- POOL_POSITIONS[base::sample(1:length(gene_pool), 
                                                       nr_random_mark, 
                                                       replace = FALSE)]
  stopifnot(iteration_df_mark[,i] %in% POOL_POSITIONS)
  stopifnot(!iteration_df_mark[,i] %in% SIGN_POSITIONS)
  stopifnot(!iteration_df_mark[,i] %in% SUBCLUSTERING_POSITIONS)
  
  colnames(iteration_df_mark)[i] <- i
}
print(iteration_df_mark[1:10,1:2])
print(dim(iteration_df_mark))

# add the signature gene positions IN THE ORIGINAL SEU DATASET
iteration_df_mark <- base::rbind(add_df_sign, iteration_df_mark)
print(iteration_df_mark[1:10,1:2])
print(dim(iteration_df_mark))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#############    all BL6 markers vs. signature + random     ####################

#-------------------------------------------------------------------------------

print("all BL6 markers vs. signature + random")

# get all BL6 marker IDs, removing subclustering genes in cell type specific way
mmus_marker_list <- lapply(geneset_list, function(geneset){
  mmus_markers <- geneset$conserved_df$gene[
    which(!is.na(geneset$conserved_df$mmus))]
  mmus_markers <- mmus_markers[
    which(!mmus_markers %in% geneset$subclustering_genes)]
  return(mmus_markers)
})
mmms_IDs <- base::unique(unlist(mmus_marker_list))

# get the number of random genes that need to be added to 
# all BL6 genes for each gene set to be permuted
nr_random_mmms <- length(mmms_IDs) - nr_sign

print(length(mmms_IDs))
print(nr_random_mmms)

#-------------------------------------------------------------------------------

set.seed(36) # different set of random numbers
base::RNGkind("L'Ecuyer-CMRG")

iteration_df_mmms <- base::data.frame(row.names = c(1:nr_random_mmms))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SCE subsetted to a smaller list 
  # can use same pool positions that 
  # - don't contain signature genes
  # - don't contain subclustering genes
  # - are expressed n >= cut_off_counts counts
  iteration_df_mmms[,i] <- POOL_POSITIONS[base::sample(1:length(gene_pool), 
                                                       nr_random_mmms, 
                                                       replace = FALSE)]
  stopifnot(iteration_df_mmms[,i] %in% POOL_POSITIONS)
  stopifnot(!iteration_df_mmms[,i] %in% SIGN_POSITIONS)
  stopifnot(!iteration_df_mmms[,i] %in% SUBCLUSTERING_POSITIONS)
  
  colnames(iteration_df_mmms)[i] <- i
}
print(iteration_df_mmms[1:10,1:2])
print(dim(iteration_df_mmms))

# add signature gene positions
iteration_df_mmms <- base::rbind(add_df_sign, iteration_df_mmms)

print(iteration_df_mmms[1:10,1:2])
print(dim(iteration_df_mmms))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

##########    all BL6 markers vs. conserved markers + random     ###############

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

##########  Prepare conserved markers IDs and seurat objects  ##################

#-------------------------------------------------------------------------------

print("all BL6 markers vs. conserved markers + random")

# conserved marker IDs have already been extracted
nr_mark <- length(mark_IDs) 
print(nr_mark)

# get position of conserved markers genes 
# POSITIONS FROM ORIGINAL SCE
MARK_POSITIONS <- which(rownames(sce_input) %in% mark_IDs)
print("MARK_POSITIONS")
print(length(MARK_POSITIONS))

print(base::table(!is.na(base::match(MARK_POSITIONS, SIGN_POSITIONS))))

#-------------------------------------------------------------------------------
# generate a fitting data frame to add conserved marker gene positions to the 
# random positions downstream

add_df_mark <- base::data.frame(row.names = c(1:nr_mark))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SCE
  add_df_mark[,i] <- MARK_POSITIONS
  colnames(add_df_mark)[i] <- i
}
print("mark_df")
print(add_df_mark[1:10,1:2])
print(dim(add_df_mark))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# GENE POOL FOR CONSERVED MARKERS

# subset a seurat object without conserved marker IDs, which contains the pool  
# of genes from which random genes can be drawn:
# - no conserved marker genes
# - no subclustering genes 
# - genes expressed at least n > cut_off_counts times in RAW counts

# check how many subclustering genes became conserved markers
print(base::table(SUBCLUSTERING_POSITIONS %in% MARK_POSITIONS))

print(length(MARK_POSITIONS))
print(length(SUBCLUSTERING_POSITIONS))
REMOVE_POS <- base::unique(c(MARK_POSITIONS, SUBCLUSTERING_POSITIONS))
print(length(REMOVE_POS))
non_mark_sce <- rownames(sce_input)[-REMOVE_POS]

sce_pool_mark <- sce_input[which(rownames(sce_input) %in% non_mark_sce),]
print(nrow(sce_input))
print(nrow(sce_pool_mark))

# get only genes that have a count of at least n = cut_off_counts in RAW counts
gene_pool_mark <- rownames(sce_pool_mark)[
  which(rowSums(assays(sce_pool_mark)$counts) >= cut_off_counts)]
print(length(gene_pool_mark))

# subset seurat object to these genes 
sce_pool_mark <- sce_pool_mark[which(rownames(sce_pool_mark) %in% 
                                       gene_pool_mark),]

print(dim(sce_pool_mark))

#-------------------------------------------------------------------------------

# the positions of genes in SEU that are allowed as pool for random genes
# POSITIONS FROM ORIGINAL SEU_PREPROCESSED
POOL_POSITIONS_MARK <- which(rownames(sce_input) %in% 
                               rownames(sce_pool_mark))

print(length(POOL_POSITIONS_MARK))

#-------------------------------------------------------------------------------

# get the number of random genes that need to be added to 
# all BL6 IDs for each gene set to be permuted for comparison with cons markers
nr_random_mmms_mark <- length(mmms_IDs) - nr_mark

print(nr_mark)
print(length(mmms_IDs))
print(nr_random_mmms_mark)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# generate i = iterations random sets of genes at the required number
# from the allowed pool of random genes (original positions in SEU)
# always generate the same random numbers
set.seed(35)
base::RNGkind("L'Ecuyer-CMRG")

# mmms_mark = comparison mmms with mark
iteration_df_mmms_mark <- base::data.frame(row.names = c(1:nr_random_mmms_mark))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SCE subsetted to a smaller list 
  # basically, subset the vector of pool positions in random positions 
  iteration_df_mmms_mark[,i] <- POOL_POSITIONS_MARK[
    base::sample(1:length(gene_pool_mark), 
                 nr_random_mmms_mark, 
                 replace = FALSE)]
  stopifnot(iteration_df_mmms_mark[,i] %in% POOL_POSITIONS_MARK)
  stopifnot(!iteration_df_mmms_mark[,i] %in% MARK_POSITIONS)
  stopifnot(!iteration_df_mmms_mark[,i] %in% SUBCLUSTERING_POSITIONS)
  
  colnames(iteration_df_mmms_mark)[i] <- i
}
print(iteration_df_mmms_mark[1:10,1:2])
print(dim(iteration_df_mmms_mark))

# add the signature gene positions IN THE ORIGINAL SEU DATASET
iteration_df_mmms_mark <- base::rbind(add_df_mark, iteration_df_mmms_mark)
print(iteration_df_mmms_mark[1:10,1:2])
print(dim(iteration_df_mmms_mark))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# TEST

# make sure there are no duplicated positions as failsave
# ALL POSITIONS FROM ORIGINAL SEU_PREPROCESSED
for(i in 1:iterations){
  stopifnot(!duplicated(iteration_df_mark[,i]))
  stopifnot(!duplicated(iteration_df_mmms[,i]))
  stopifnot(!duplicated(iteration_df_mmms_mark[,i]))
}

# stichproben
# signature genes should be in all 

# in two columns of the same iteration DF
print(base::paste("should be at least", nr_sign, "TRUE"))
print(base::table(!is.na(base::match(iteration_df_mark[,1],
                                     iteration_df_mark[,2]))))

# in two columns of the mark and mmms iteration DFs
print(base::paste("should be at least", nr_sign, "TRUE"))
print(base::table(!is.na(base::match(iteration_df_mark[,1],
                                     iteration_df_mmms[,2]))))

# conserved markers should be in two columsn of the same iteration DF
print(base::paste("should be at least", nr_mark, "TRUE"))
print(base::table(!is.na(base::match(iteration_df_mmms_mark[,1],
                                     iteration_df_mmms_mark[,2]))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# run standard seurat re-clustering pipeline (own function) i times

print("starting iteration mark vs. signature + random")

res_df_list_mark <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = permuting_reclustering_scores_sce,
  sce = sce_input, # use original sce so positions are correct
  iteration_df = iteration_df_mark,
  resolution_louvain = resolution_louvain,
  k_graph = k_graph,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

# res_df_list_mark <- lapply(
#   X = as.list(c(1:iterations)),
#   FUN = permuting_reclustering_scores_sce,
#   sce = sce_input, # use original sce
#   resolution_louvain = resolution_louvain,
#   k_graph = k_graph,
#   iteration_df = iteration_df_mark)

score_df_mark <- dplyr::bind_rows(res_df_list_mark)
print(head(score_df_mark))

print("starting iteration mmms vs. signature + random")

res_df_list_mmms <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = permuting_reclustering_scores_sce,
  sce = sce_input,
  resolution_louvain = resolution_louvain,
  k_graph = k_graph,
  iteration_df = iteration_df_mmms,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

score_df_mmms <- dplyr::bind_rows(res_df_list_mmms)
print(head(score_df_mmms))

print("starting iteration mmms vs. conserved markers + random")

res_df_list_mmms_mark <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = permuting_reclustering_scores_sce,
  sce = sce_input,
  resolution_louvain = resolution_louvain,
  k_graph = k_graph,
  iteration_df = iteration_df_mmms_mark,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE)

score_df_mmms_mark <- dplyr::bind_rows(res_df_list_mmms_mark)
print(head(score_df_mmms_mark))

#-------------------------------------------------------------------------------

base::saveRDS(score_df_mark, snakemake@output[["perm_score_df_mark"]])
base::saveRDS(score_df_mmms, snakemake@output[["perm_score_df_mmms"]])

base::saveRDS(score_df_mmms_mark, snakemake@output[["perm_score_df_mmms_mark"]])

utils::sessionInfo()
