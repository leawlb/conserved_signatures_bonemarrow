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

# params
cons_level_use <- snakemake@params[["cons_level_use"]]
comparison <- snakemake@params[["comparison"]]
dataset_curr <- snakemake@wildcards[["dataset"]]

# alternative shorter word used as name for the list items
cons_level_use_alt <- snakemake@params[["cons_level_use_alt"]]

print(comparison)
print(cons_level_use)
print(dataset_curr)

#-------------------------------------------------------------------------------

# list with dfs of original scores to be tested, for each conservation level
orig_score_df_list <- base::readRDS(snakemake@input[["orig_score_df_list_input"]])

# df of permutated scores from the conservation level and comp as indicated
perm_score_df <- base::readRDS(snakemake@input[["perm_score_df_input"]])

print(head(perm_score_df))

#-------------------------------------------------------------------------------

# resolution DF
# load the dataframe which contains info on which resolution to use for
# the conservation level as indicated
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
resl <- resolution_df$resolution[
  resolution_df$conservation_level == cons_level_use]

print(resl)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# subset to correct orig data frame
orig_score_df <- orig_score_df_list[[cons_level_use_alt]][[resl]]

# add info and subset
orig_score_df$dataset <- base::rep(dataset_curr, nrow(orig_score_df))
orig_score_df$comparison <- base::rep(comparison, nrow(orig_score_df))

orig_score_df$nr_genes_used <- base::rep(perm_score_df$nr_genes_used[1], 
                                         nrow(orig_score_df))

colnames(orig_score_df)[colnames(orig_score_df) == "value"] <- "value_orig"

# failsave
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
# correct pvals
# using bonferroni hochberg for now
# possibility of merging all pval_score_df_output and correcting pvals together

# orig_score_df_sub$pval_corrected <- stats::p.adjust(
#   orig_score_df_sub$pval, 
#   method = "BH")

#-------------------------------------------------------------------------------

base::saveRDS(orig_score_df_sub, snakemake@output[["pval_score_df_output"]])
