#-------------------------------------------------------------------------------

# export pvals from permutations and correct them

# this script will be run for multiple comparisons:
# background:
# - conserved signature vs random ("sign-vs-rand")

# geneset comparison:
# - conserved markers vs signature + random ("mark-vs-signrand")
# - all BL6 markers vs conserved signature + random ("mmms-vs-signrand")
# - all BL6 markers vs conserved markers + random ("mmms-vs-markrand")

#-------------------------------------------------------------------------------

# determine random number generator for sample
# "Mersenne-Twister" is default
RNGkind("Mersenne-Twister") 
set.seed(37)

#-------------------------------------------------------------------------------
# load

# df of original scores to be tested
orig_score_df <- base::readRDS(snakemake@input[["orig_score_df_input"]])
# remove fraction columns for downstream compatibility with other datasets
orig_score_df <- orig_score_df[,-which(colnames(orig_score_df) == "fraction")]

# df of permutated scores from the same object
perm_score_df <- base::readRDS(snakemake@input[["perm_score_df_input"]])

# params
cons_level_use <- snakemake@params[["cons_level_use"]]
comparison <- snakemake@params[["comparison"]]
fraction_curr <- snakemake@wildcards[["fraction"]]

# get info on the resolution to record it in the DF
resolution_louvain_list <- snakemake@params[["resolution_louvain_list"]]
resolution <- resolution_louvain_list[[fraction_curr]]

# add info and subset
# call condition so it fits with 02
print(head(orig_score_df))
orig_score_df$condition <- base::rep(fraction_curr, nrow(orig_score_df))
orig_score_df$comparison <- base::rep(comparison, nrow(orig_score_df))
orig_score_df$identifier <- base::rep("own", nrow(orig_score_df))

orig_score_df_sub <- orig_score_df[
  orig_score_df$conservation_level %in% cons_level_use,]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# calculate P val based on original score and permutated scores for each 
# type of score

#-------------------------------------------------------------------------------
# list of scores to iterate through

scores_v <- c(
  "adjusted_rand_index",
  "mean_prop_cells_cluster",
  "variation_information")

#-------------------------------------------------------------------------------
# per score

orig_score_df_sub <- orig_score_df_sub[orig_score_df_sub$type %in% scores_v,]
orig_score_df_sub$nr_iterations <- vector(length = nrow(orig_score_df_sub))
orig_score_df_sub$median_perm_scores <- vector(length = nrow(orig_score_df_sub))
orig_score_df_sub$pval <- vector(length = nrow(orig_score_df_sub))

for(score in scores_v){
  
  print(score)
  
  # further subset orig_score_df_sub
  orig_score_df_temp <- orig_score_df_sub[orig_score_df_sub$type == score,]
  perm_scores_all <- perm_score_df[,which(colnames(perm_score_df) == score)]
  
  # the pval is the proportion of permuted values that is at least as "good"
  # as the original value 
  if(score != "variation_information"){
    pval = length(which(perm_scores_all >= orig_score_df_temp$value))/
      length(perm_scores_all)
  }else if(score == "variation_information"){
    # for variation_information score, the smaller, the better
    pval = length(which(perm_scores_all <= orig_score_df_temp$value))/
      length(perm_scores_all)
  }
  
  median_perm_scores <- stats::median(perm_scores_all)
  
  orig_score_df_sub$pval[
    orig_score_df_sub$type == score] <- pval
  orig_score_df_sub$median_perm_scores[
    orig_score_df_sub$type == score] < median_perm_scores
  orig_score_df_sub$nr_iterations[
    orig_score_df_sub$type == score] <- length(perm_scores_all)
  
}

print(head(orig_score_df_sub))

#-------------------------------------------------------------------------------
# add resolution

print(base::paste("Resolution", resolution))
orig_score_df_sub$resolution <- base::rep(resolution, nrow(orig_score_df_sub))

#-------------------------------------------------------------------------------
# change name of value column

colnames(orig_score_df_sub)[colnames(orig_score_df_sub) == "value"] <- "value_orig"
print(colnames(orig_score_df_sub))

# make sure the fraction column is gone
stopifnot(!"fraction" %in% colnames(orig_score_df_sub))
stopifnot("resolution" %in% colnames(orig_score_df_sub))

#-------------------------------------------------------------------------------

base::saveRDS(orig_score_df_sub, snakemake@output[["pval_score_df_output"]])
