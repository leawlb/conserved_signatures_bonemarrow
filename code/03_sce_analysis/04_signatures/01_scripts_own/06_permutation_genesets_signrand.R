
# permutation to compare genes sets (conserved markers and all BL6 markers)
# vs. conserved signature + random genes (same nr) to see if they contain
# important genes in addition to containing conserved signature genes.

# mmms_mark = BL6-only markers vs conserved markers + random
# these steps have been # commented because the comparison is not needed
# but the code is still kept.
# the respective lines from the snakefile have been deleted, though

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# determine random number generator for sample
library(parallel, quietly = TRUE)
RNGkind("L'Ecuyer-CMRG") # using this for parallel is necessary
set.seed(37)

#-------------------------------------------------------------------------------

library(dplyr, quietly = TRUE)
library(BiocGenerics, quietly = TRUE)
library(S4Vectors, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)

# for calculating scores
library(mclust, quietly = TRUE)
library(bluster, quietly = TRUE)

source(snakemake@params[["functions_reclustering"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD OBJECTS

# params for re-clustering and permutation test
cut_off_prop <- snakemake@params[["cut_off_prop"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

#-------------------------------------------------------------------------------

# load the seurat object/dataset to be tested
# - subclustering genes removed
# - cts to exclude removed
# - hsc dataset sub-sampled
sce_input <- base::readRDS(snakemake@input[["sce_input"]])
print(dim(sce_input))

#-------------------------------------------------------------------------------
# params for re-clustering

k_graph_list <- snakemake@params[["k_graph_list"]]
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
fraction_curr <-  snakemake@wildcards[["fraction"]]

k_graph <- k_graph_list[[fraction_curr]]
resolution_louvain <- resolution_louvain_list[[fraction_curr]]

#-------------------------------------------------------------------------------
# load geneset list
geneset_list <- base::readRDS(snakemake@input[["geneset_list"]])

#-------------------------------------------------------------------------------

# set seeds for sample()

if(fraction_curr == "hsc"){
  seed1 <- 33
  seed2 <- 32
  seed3 <- 31
}else if(fraction_curr == "str"){
  seed1 <- 999
  seed2 <- 9999
  seed3 <- 99999
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

################   Prepare signature IDs and sce object    #####################

#-------------------------------------------------------------------------------

print("preparing signature IDs")

# get all signature genes, removing subclustering genes in ct specific way
conserved_signature_list <- lapply(geneset_list, function(geneset){
  conserved_signature <- geneset$conserved_signature
  return(conserved_signature)
})
sign_IDs <- base::unique(unlist(conserved_signature_list))
sign_IDs <- sign_IDs[sign_IDs %in% rownames(sce_input)]

nr_sign <- length(sign_IDs) 
print("nr_sign")
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
print(dim(add_df_sign))

#-------------------------------------------------------------------------------
# extract subclustering genes and control

subclustering_gene_list <- lapply(geneset_list, function(geneset){
  subclustering_genes <- geneset$genes_subclustering
  return(subclustering_genes)
})
subclustering_genes <- base::unique(unlist(subclustering_gene_list))
stopifnot(!subclustering_genes %in% rownames(sce_input))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

############   Conserved markers vs. signature + random     ####################

#-------------------------------------------------------------------------------

print("Conserved markers vs. signature + random")

# extract conserved marker IDs
# keep only genes that are also in SCE object, which already has 
# sub-clustering genes removed in total
conserved_markers_list <- lapply(geneset_list, function(geneset){
  conserved_markers <- geneset$conserved_markers
  return(conserved_markers)
})
mark_IDs <- base::unique(unlist(conserved_markers_list))
mark_IDs <- mark_IDs[mark_IDs %in% rownames(sce_input)]

# get the number of random genes that need to be added to 
# conserved_signature_IDs for each gene set to be permuted
nr_random_mark <- length(mark_IDs) - nr_sign

print("mark_IDs")
print(length(mark_IDs))
print(nr_random_mark)

#-------------------------------------------------------------------------------

# subset a SCE object without signature IDs, which contains the pool of 
# genes from which random genes can be drawn:
# - no signature genes
# - no subclustering genes
# - genes expressed at least >= cut_off_prop of cells (in format 0.01 = 1%) 

print("subset")
print(nrow(sce_input))
print(length(SIGN_POSITIONS))
non_sign_sce <- rownames(sce_input)[-SIGN_POSITIONS]

# subset SCE object to obtain non-conserved signature genes
sce_pool <- sce_input[which(rownames(sce_input) %in% non_sign_sce),]

# get only genes that are expressed in at least cut_off_prop% of cells
# use own function to get the proportion of cells a gene is expressed in
gene_pool <- rownames(sce_pool)
print(length(gene_pool))

prop_df <- prop_expressed_total_sce(
  sce = sce_pool, 
  geneset = gene_pool)

# subset by cut-off proportion
prop_df_sub <- prop_df[prop_df$prop_cells >= cut_off_prop,]

gene_pool <- gene_pool[which(gene_pool %in% prop_df_sub$gene)]
print("after cut-off")
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

set.seed(seed1) # different set of random numbers each time
for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SCE subsetted to a smaller list 
  # basically, subset the vector of pool positions in random positions 
  iteration_df_mark[,i] <- POOL_POSITIONS[base::sample(1:length(gene_pool), 
                                                       nr_random_mark, 
                                                       replace = FALSE)]
  stopifnot(iteration_df_mark[,i] %in% POOL_POSITIONS)
  stopifnot(!iteration_df_mark[,i] %in% SIGN_POSITIONS)

  colnames(iteration_df_mark)[i] <- i
}
set.seed(37)

print(iteration_df_mark[1:10,1:2])

# add the signature gene positions IN THE ORIGINAL SEU DATASET
iteration_df_mark <- base::rbind(add_df_sign, iteration_df_mark)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#############    all BL6 markers vs. signature + random     ####################

#-------------------------------------------------------------------------------

print("all BL6 markers vs. signature + random")

# get all BL6 marker IDs
# keep only genes that are also in SCE object, which already has 
# sub-clustering genes removed in total
mmus_marker_list <- lapply(geneset_list, function(geneset){
  mmus_markers <- geneset$conserved_df$gene[
    which(!is.na(geneset$conserved_df$mmus))]
  return(mmus_markers)
})
mmms_IDs <- base::unique(unlist(mmus_marker_list))
mmms_IDs <- mmms_IDs[mmms_IDs %in% rownames(sce_input)]

# get the number of random genes that need to be added to 
# all BL6 genes for each gene set to be permuted
nr_random_mmms <- length(mmms_IDs) - nr_sign

print("mmms_IDs")
print(length(mmms_IDs))
print(nr_random_mmms)

#-------------------------------------------------------------------------------

iteration_df_mmms <- base::data.frame(row.names = c(1:nr_random_mmms))

set.seed(seed2) 
for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SCE subsetted to a smaller list 
  # can use same pool positions that 
  # - don't contain signature genes
  # - don't contain subclustering genes
  # - are expressed n >= cut_off_prop% of cells
  iteration_df_mmms[,i] <- POOL_POSITIONS[base::sample(1:length(gene_pool), 
                                                       nr_random_mmms, 
                                                       replace = FALSE)]
  stopifnot(iteration_df_mmms[,i] %in% POOL_POSITIONS)
  stopifnot(!iteration_df_mmms[,i] %in% SIGN_POSITIONS)

  colnames(iteration_df_mmms)[i] <- i
}
set.seed(37) 

print(iteration_df_mmms[1:10,1:2])

# add signature gene positions
iteration_df_mmms <- base::rbind(add_df_sign, iteration_df_mmms)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

##########    all BL6 markers vs. conserved markers + random     ###############

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

##########  Prepare conserved markers IDs and seurat objects  ##################

#-------------------------------------------------------------------------------

# print("all BL6 markers vs. conserved markers + random")
# 
# # conserved marker IDs have already been extracted
# nr_mark <- length(mark_IDs) 
# 
# # get position of conserved markers genes 
# # POSITIONS FROM ORIGINAL SCE
# MARK_POSITIONS <- which(rownames(sce_input) %in% mark_IDs)
# print("MARK_POSITIONS")
# print(length(MARK_POSITIONS))
# print(base::table(!is.na(base::match(MARK_POSITIONS, SIGN_POSITIONS))))

#-------------------------------------------------------------------------------
# generate a fitting data frame to add conserved marker gene positions to the 
# random positions downstream

# add_df_mark <- base::data.frame(row.names = c(1:nr_mark))
# 
# for(i in 1:iterations){
#   # POSITIONS FROM ORIGINAL SCE
#   add_df_mark[,i] <- MARK_POSITIONS
#   colnames(add_df_mark)[i] <- i
# }

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# GENE POOL FOR CONSERVED MARKERS

# subset a seurat object without conserved marker IDs, which contains the pool  
# of genes from which random genes can be drawn:
# - no subclustering genes (already removed)
# - no conserved marker genes
# - genes expressed expressed in >= cut_off_prop % of cells

# remove conserved marker genes
# non_mark_sce <- rownames(sce_input)[-MARK_POSITIONS]
# 
# sce_pool_mark <- sce_input[which(rownames(sce_input) %in% non_mark_sce),]
# gene_pool_mark <- rownames(sce_pool_mark)
# print("subset")
# print(nrow(sce_input))
# print(length(gene_pool_mark))
# 
# # remove genes expressed in < cut_off_prop % of cells
# prop_df_mark <- prop_expressed_total_sce(
#   sce = sce_pool_mark, 
#   geneset = gene_pool_mark)
# 
# # subset by cut-off proportion
# prop_df_sub_mark <- prop_df_mark[prop_df_mark$prop_cells >= cut_off_prop,]
# print("after cut-off")
# print(nrow(prop_df_sub_mark))
# 
# gene_pool_mark <- gene_pool_mark[which(gene_pool_mark %in%
#                                          prop_df_sub_mark$gene)]
# 
# # subset seurat object to these genes 
# sce_pool_mark <- sce_pool_mark[which(rownames(sce_pool_mark) %in% 
#                                        gene_pool_mark),]
# 
# print(dim(sce_pool_mark))

#-------------------------------------------------------------------------------

# # the positions of genes in SEU that are allowed as pool for random genes
# # POSITIONS FROM ORIGINAL SEU_PREPROCESSED
# POOL_POSITIONS_MARK <- which(rownames(sce_input) %in% 
#                                rownames(sce_pool_mark))
# 
# print(length(POOL_POSITIONS_MARK))

#-------------------------------------------------------------------------------

# # get the number of random genes that need to be added to 
# # all BL6 IDs for each gene set to be permuted for comparison with cons markers
# # nr_random_mmms_mark <- length(mmms_IDs) - nr_mark
# 
# print("nr_random_mmms_mark")
# print(length(mmms_IDs))
# print(nr_mark)
# print(nr_random_mmms_mark)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# # generate i = iterations random sets of genes at the required number
# # from the allowed pool of random genes (original positions in SEU)
# # always generate the same random numbers
# 
# # mmms_mark = comparison mmms with mark
# iteration_df_mmms_mark <- base::data.frame(row.names = c(1:nr_random_mmms_mark))
# 
# set.seed(seed3)
# for(i in 1:iterations){
#   iteration_df_mmms_mark[,i] <- POOL_POSITIONS_MARK[
#     base::sample(1:length(gene_pool_mark),
#                  nr_random_mmms_mark,
#                  replace = FALSE)]
#   stopifnot(iteration_df_mmms_mark[,i] %in% POOL_POSITIONS_MARK)
#   stopifnot(!iteration_df_mmms_mark[,i] %in% MARK_POSITIONS)
# 
#   colnames(iteration_df_mmms_mark)[i] <- i
# }
# set.seed(37)
# 
# print(iteration_df_mmms_mark[1:10, 1:2])
# 
# print("iteration_df_mmms_mark before addition of SIGN_POS")
# print(dim(iteration_df_mmms_mark))
# # add the signature gene positions IN THE ORIGINAL SEU DATASET
# iteration_df_mmms_mark <- base::rbind(add_df_mark, iteration_df_mmms_mark)
# print("iteration_df_mmms_mark after addition of SIGN_POS")
# print(dim(iteration_df_mmms_mark))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# TEST

# make sure there are no duplicated positions as failsave
# ALL POSITIONS FROM ORIGINAL SEU_PREPROCESSED
for(i in 1:iterations){
  stopifnot(!duplicated(iteration_df_mark[,i]))
  stopifnot(!duplicated(iteration_df_mmms[,i]))
  #stopifnot(!duplicated(iteration_df_mmms_mark[,i]))
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

# # conserved markers should be in two columsn of the same iteration DF
# print(base::paste("should be at least", nr_mark, "TRUE"))
# print(base::table(!is.na(base::match(iteration_df_mmms_mark[,1],
#                                      iteration_df_mmms_mark[,2]))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ITERATION

print("nr_cores")
print(nr_cores)
print("iterations")
print(iterations)
print("cut_off_prop")
print(cut_off_prop)
print("k_graph")
print(k_graph)
print("resolution_louvain")
print(resolution_louvain)

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
  mc.silent = TRUE,
  mc.set.seed = TRUE)

score_df_mark <- dplyr::bind_rows(res_df_list_mark)
score_df_mark$cut_off_prop <- base::rep(cut_off_prop, nrow(score_df_mark))
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
  mc.silent = TRUE,
  mc.set.seed = TRUE)

score_df_mmms <- dplyr::bind_rows(res_df_list_mmms)
score_df_mmms$cut_off_prop <- base::rep(cut_off_prop, nrow(score_df_mmms))
print(head(score_df_mmms))

# print("starting iteration mmms vs. conserved markers + random")
# 
# res_df_list_mmms_mark <- parallel::mclapply(
#   X = as.list(c(1:iterations)),
#   FUN = permuting_reclustering_scores_sce,
#   sce = sce_input,
#   resolution_louvain = resolution_louvain,
#   k_graph = k_graph,
#   iteration_df = iteration_df_mmms_mark,
#   mc.preschedule = TRUE,
#   mc.cores = nr_cores,
#   mc.silent = TRUE,
#   mc.set.seed = TRUE)
# 
# score_df_mmms_mark <- dplyr::bind_rows(res_df_list_mmms_mark)
# score_df_mmms_mark$cut_off_prop <- base::rep(cut_off_prop, 
#                                              nrow(score_df_mmms_mark))
# print(head(score_df_mmms_mark))

#-------------------------------------------------------------------------------

# for testing without mclapply
# res_df_list_mark <- lapply(
#   X = as.list(c(1:iterations)),
#   FUN = permuting_reclustering_scores_sce,
#   sce = sce_input, # use original sce
#   resolution_louvain = resolution_louvain,
#   k_graph = k_graph,
#   iteration_df = iteration_df_mark)

# score_df_mark <- dplyr::bind_rows(res_df_list_mark)
# score_df_mark$cut_off_prop <- base::rep(cut_off_prop, nrow(score_df_mark))
# print(head(score_df_mark))

#-------------------------------------------------------------------------------

base::saveRDS(score_df_mark, snakemake@output[["perm_score_df_mark"]])
base::saveRDS(score_df_mmms, snakemake@output[["perm_score_df_mmms"]])

# this object is not exported, because it is not needed
#base::saveRDS(score_df_mmms_mark, snakemake@output[["perm_score_df_mmms_mark"]])

utils::sessionInfo()
