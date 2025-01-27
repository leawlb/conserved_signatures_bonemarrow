
# EXPLANATION TODO


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# determine random number generator for sample()
# determine random number generator for sample
library(parallel)
RNGkind("L'Ecuyer-CMRG") # using this for usages of parallel is necessary

#-------------------------------------------------------------------------------

library(Seurat, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(BiocGenerics, quietly = TRUE)
library(S4Vectors, quietly = TRUE)

# for calculating scores
library(mclust, quietly = TRUE)
library(bluster, quietly = TRUE)

source(snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# LOAD OBJECTS

# params for re-clustering and permutation test
cut_off_prop <- snakemake@params[["cut_off_prop"]]
nr_cores <- snakemake@params[["nr_cores"]]
iterations <- snakemake@params[["iterations"]]

#-------------------------------------------------------------------------------

# determine which gene sets to use for the current dataset (based on fraction)
dataset_curr <- snakemake@wildcards[["dataset"]]

# determine the current fraction to choose the correct gene sets
datasets_other_hsc <- snakemake@params[["datasets_other_hsc"]]
datasets_other_str <- snakemake@params[["datasets_other_str"]]

if(dataset_curr %in% datasets_other_str){
  fraction_curr <- "str"
}else if(dataset_curr %in% datasets_other_hsc){
  fraction_curr <- "hsc"
}
print(fraction_curr)

#-------------------------------------------------------------------------------

# load the seurat object/dataset to be tested
seu_preprocessed <- base::readRDS(snakemake@input[["seu_preprocessed"]])

#-------------------------------------------------------------------------------

# load ensembl IDs of specific gene set to be tested from correct fraction
ensembl_sign <- snakemake@input[["ensembl_sign"]]
ensembl_mark <- snakemake@input[["ensembl_mark"]]
ensembl_mmms <- snakemake@input[["ensembl_mmms"]]

ensembl_sign_df <- base::readRDS(ensembl_sign[[
  which(base::grepl(fraction_curr, ensembl_sign))]])
ensembl_mark_df <- base::readRDS(ensembl_mark[[
  which(base::grepl(fraction_curr, ensembl_mark))]])
ensembl_mmms_df <- base::readRDS(ensembl_mmms[[
  which(base::grepl(fraction_curr, ensembl_mmms))]])

#  which gene IDs to use
ensembl_column_use <- seu_preprocessed@misc$ensembl_column_use
print(dataset_curr)
print("ensembl_column_use")
print(ensembl_column_use)

#-------------------------------------------------------------------------------

# load the dataframe which contains info on which resolution to use for
resolution_df_path <- snakemake@params[["resolution_df"]] 

resolution_df <- utils::read.csv(file = resolution_df_path, 
                                 header = TRUE, 
                                 sep = ";", 
                                 check.names=FALSE, 
                                 stringsAsFactors=FALSE, 
                                 as.is=TRUE, 
                                 colClasses = "character")

resolution_df <- resolution_df[resolution_df$dataset == dataset_curr,]

# get resolution for each gene set to be permuted
resl_mark <- resolution_df$resolution[
  resolution_df$conservation_level == "conserved_markers"]

resl_mmms <- resolution_df$resolution[
  resolution_df$conservation_level == "mmusall_markers"]

print("resolutions:")
print(resl_mark)
print(resl_mmms)

#-------------------------------------------------------------------------------
# determine seeds

seed1 <- (99 + round(ncol(seu_preprocessed)/100, digits = 0))
seed2 <- (98 + round(ncol(seu_preprocessed)/100, digits = 0))
seed3 <- (97 + round(ncol(seu_preprocessed)/100, digits = 0))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

############   Prepare signature IDs and seurat objects     ####################

#-------------------------------------------------------------------------------

print("preparing signature IDs")

# extract conserved signature IDs, keep only genes also in seu object
sign_IDs <- base::unique(
  ensembl_sign_df[,which(colnames(ensembl_sign_df) == ensembl_column_use)])
sign_IDs <- sign_IDs[sign_IDs %in% rownames(seu_preprocessed)]

nr_sign <- length(sign_IDs) 
print("nr_sign")
print(nr_sign)

#-------------------------------------------------------------------------------
# get position of conserved signature genes 
# POSITIONS FROM ORIGINAL SEU_PREPROCESSED

SIGN_POSITIONS <- which(rownames(seu_preprocessed) %in% sign_IDs)
print(length(SIGN_POSITIONS))

#-------------------------------------------------------------------------------
# generate a fitting data frame to add signature gene positions to the 
# random positions downstream

add_df_sign <- base::data.frame(row.names = c(1:nr_sign))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED
  add_df_sign[,i] <- SIGN_POSITIONS
  colnames(add_df_sign)[i] <- i
}
print("ensembl_sign_df")
print(add_df_sign[1:10,1:2])
print(nrow(add_df_sign))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

############   Conserved markers vs. signature + random     ####################

#-------------------------------------------------------------------------------

print("Conserved markers vs. signature + random")

# extract conserved marker IDs
mark_IDs <- base::unique(
  ensembl_mark_df[,which(colnames(ensembl_mark_df) == ensembl_column_use)])
mark_IDs <- mark_IDs[mark_IDs %in% rownames(seu_preprocessed)]

# get the number of random genes that need to be added to 
# conserved_signature_IDs for each gene set to be permuted
nr_random_mark <- length(mark_IDs) - nr_sign

print(length(mark_IDs))
print(nr_sign)
print(nr_random_mark)
print(nr_random_mark + nr_sign)

#-------------------------------------------------------------------------------

# subset a seurat object without signature IDs, which contains the pool of 
# genes from which random genes can be drawn:
# - no signature genes
# - genes expressed in at least cut_off_prop cells 

non_seu_sign <- rownames(seu_preprocessed)[-SIGN_POSITIONS]

seu_pool <- BiocGenerics::subset(seu_preprocessed,
                                 features = non_seu_sign,
                                 slot = "count")

# get only genes that are expressed in at least cut_off_prop% of cells
# use own function to get the proportion of cells a gene is expressed in
gene_pool <- rownames(seu_pool)
print(length(gene_pool))

prop_df <- prop_expressed_total_seu(
  seu = seu_pool, 
  geneset = gene_pool)

# subset by cut-off proportion
prop_df_sub <- prop_df[prop_df$prop_cells >= cut_off_prop,]

gene_pool <- gene_pool[which(gene_pool %in% prop_df_sub$gene)]
print(length(gene_pool))

# subset seurat object to these genes 
seu_pool <- BiocGenerics::subset(seu_pool,
                                 features = gene_pool,
                                 slot = "count")
print(dim(seu_pool))
  
#-------------------------------------------------------------------------------

# the positions of genes in SEU that are allowed as pool for random genes
# POSITIONS FROM ORIGINAL SEU_PREPROCESSED
POOL_POSITIONS <- which(rownames(seu_preprocessed) %in% rownames(seu_pool))
print(length(POOL_POSITIONS))

#-------------------------------------------------------------------------------

# generate i = iterations random sets of genes at the required number
# from the allowed pool of random genes (original positions in SEU)
# always generate the same random numbers

iteration_df_mark <- base::data.frame(row.names = c(1:nr_random_mark))

# make sure that seed is unique for each object
set.seed(seed1)
for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED subsetted to a smaller list 
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
print(dim(iteration_df_mark))

# add the signature gene positions IN THE ORIGINAL SEU DATASET
iteration_df_mark <- base::rbind(add_df_sign, iteration_df_mark)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#############    all BL6 markers vs. signature + random     ####################

#-------------------------------------------------------------------------------

print("all BL6 markers vs. signature + random")

# get all BL6 marker IDs
mmms_IDs <- base::unique(
  ensembl_mmms_df[,which(colnames(ensembl_mmms_df) == ensembl_column_use)])
mmms_IDs <- mmms_IDs[mmms_IDs %in% rownames(seu_preprocessed)]

# get the number of random genes that need to be added to 
# all BL6 genes for each gene set to be permuted
nr_random_mmms <- length(mmms_IDs) - nr_sign

print(length(mmms_IDs))
print(nr_random_mmms)

#-------------------------------------------------------------------------------


iteration_df_mmms <- base::data.frame(row.names = c(1:nr_random_mmms))

set.seed(seed2)
for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED subsetted to a smaller list 
  # can use same pool positions that 
  # - don't contain signature genes
  # - are expressed in at least cut_off_prop % of cells
  iteration_df_mmms[,i] <- POOL_POSITIONS[base::sample(1:length(gene_pool), 
                                                       nr_random_mmms, 
                                                       replace = FALSE)]
  stopifnot(iteration_df_mmms[,i] %in% POOL_POSITIONS)
  stopifnot(!iteration_df_mmms[,i] %in% SIGN_POSITIONS)
  
  colnames(iteration_df_mmms)[i] <- i
}
set.seed(37)

print(iteration_df_mmms[1:10,1:2])
print(dim(iteration_df_mmms))

# add signature gene positions
iteration_df_mmms <- base::rbind(add_df_sign, iteration_df_mmms)

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
# POSITIONS FROM ORIGINAL SEU_PREPROCESSED
MARK_POSITIONS <- which(rownames(seu_preprocessed) %in% mark_IDs)
print("MARK_POSITIONS")
print(length(MARK_POSITIONS))

print(base::table(!is.na(base::match(MARK_POSITIONS, SIGN_POSITIONS))))

#-------------------------------------------------------------------------------
# generate a fitting data frame to add conserved marker gene positions to the 
# random positions downstream

add_df_mark <- base::data.frame(row.names = c(1:nr_mark))

for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED
  add_df_mark[,i] <- MARK_POSITIONS
  colnames(add_df_mark)[i] <- i
}
print("mark_df")
print(dim(add_df_mark))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# GENE POOL FOR CONSERVED MARKERS

# subset a seurat object without conserved marker IDs, which contains the pool  
# of genes from which random genes can be drawn:
# - no conserved marker genes
# - genes expressed in at least cut_off_prop cells 

non_seu_mark <- rownames(seu_preprocessed)[-MARK_POSITIONS]

seu_pool_mark <- BiocGenerics::subset(seu_preprocessed,
                                      features = non_seu_mark,
                                      slot = "count")

# get only genes that are expressed in at least cut_off_prop% of cells
# use own function to get the proportion of cells a gene is expressed in
gene_pool_mark <- rownames(seu_pool_mark)
print(length(gene_pool_mark))

prop_df <- prop_expressed_total_seu(
  seu = seu_pool_mark, 
  geneset = gene_pool_mark)

# subset by cut-off proportion
prop_df_sub <- prop_df[prop_df$prop_cells >= cut_off_prop,]

gene_pool_mark <- gene_pool_mark[which(gene_pool_mark %in% prop_df_sub$gene)]
print(length(gene_pool_mark))

# subset seurat object to these genes 
seu_pool_mark <- BiocGenerics::subset(seu_pool_mark,
                                      features = gene_pool_mark,
                                      slot = "count")
print(dim(seu_pool_mark))

#-------------------------------------------------------------------------------

# the positions of genes in SEU that are allowed as pool for random genes
# POSITIONS FROM ORIGINAL SEU_PREPROCESSED
POOL_POSITIONS_MARK <- which(rownames(seu_preprocessed) %in% 
                               rownames(seu_pool_mark))

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

# mmms_mark = comparison mmms with mark
iteration_df_mmms_mark <- base::data.frame(row.names = c(1:nr_random_mmms_mark))

set.seed(seed3)
for(i in 1:iterations){
  # POSITIONS FROM ORIGINAL SEU_PREPROCESSED subsetted to a smaller list 
  # basically, subset the vector of pool positions in random positions 
  iteration_df_mmms_mark[,i] <- POOL_POSITIONS_MARK[
    base::sample(1:length(gene_pool_mark), 
                 nr_random_mmms_mark, 
                 replace = FALSE)]
  stopifnot(iteration_df_mmms_mark[,i] %in% POOL_POSITIONS_MARK)
  stopifnot(!iteration_df_mmms_mark[,i] %in% MARK_POSITIONS)
  
  colnames(iteration_df_mmms_mark)[i] <- i
}
set.seed(37)

print(iteration_df_mmms_mark[1:10,1:2])
print(dim(iteration_df_mmms_mark))

# add the signature gene positions IN THE ORIGINAL SEU DATASET
iteration_df_mmms_mark <- base::rbind(add_df_mark, iteration_df_mmms_mark)

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
  FUN = permuting_reclustering_scores_seurat,
  seu = seu_preprocessed,
  data_use = seu_preprocessed@misc$data_use,
  iteration_df = iteration_df_mark,
  resolution = resl_mark,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE)

# for testing without mclapply
# res_df_list_mark <- lapply(
#   X = as.list(c(1:iterations)),
#   FUN = permuting_reclustering_scores_seurat,
#   seu = seu_preprocessed,
#   assay_use = "RNA",
#   iteration_df = iteration_df_mark,
#   resolution = resl_mark)

score_df_mark <- dplyr::bind_rows(res_df_list_mark)
print(head(score_df_mark))

print("starting iteration mmms vs. signature + random")

res_df_list_mmms <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = permuting_reclustering_scores_seurat,
  seu = seu_preprocessed,
  data_use = seu_preprocessed@misc$data_use,
  iteration_df = iteration_df_mmms,
  resolution = resl_mmms,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE)

score_df_mmms <- dplyr::bind_rows(res_df_list_mmms)
print(head(score_df_mmms))

print("starting iteration mmms vs. conserved markers + random")

res_df_list_mmms_mark <- parallel::mclapply(
  X = as.list(c(1:iterations)),
  FUN = permuting_reclustering_scores_seurat,
  seu = seu_preprocessed,
  data_use = seu_preprocessed@misc$data_use,
  iteration_df = iteration_df_mmms_mark,
  resolution = resl_mmms,
  mc.preschedule = TRUE,
  mc.cores = nr_cores,
  mc.silent = TRUE,
  mc.set.seed = TRUE)

score_df_mmms_mark <- dplyr::bind_rows(res_df_list_mmms_mark)
print(head(score_df_mmms_mark))

#-------------------------------------------------------------------------------

base::saveRDS(score_df_mark, snakemake@output[["perm_score_df_mark"]])
base::saveRDS(score_df_mmms, snakemake@output[["perm_score_df_mmms"]])
base::saveRDS(score_df_mmms_mark, snakemake@output[["perm_score_df_mmms_mark"]])

utils::sessionInfo()
