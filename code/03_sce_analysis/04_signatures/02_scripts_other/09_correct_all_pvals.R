
#-------------------------------------------------------------------------------

# import ALL computed pvals and correct them

RNGkind("L'Ecuyer-CMRG")
set.seed(37)

library(dplyr)

#-------------------------------------------------------------------------------
# get all paths into one list for each main branch (own vs other)

list_own <- unlist(list(
  snakemake@input[["own_sign_rand"]], 
  snakemake@input[["own_mark_rand"]], 
  snakemake@input[["own_mmms_rand"]], 
  snakemake@input[["own_mark_signrand"]],
  snakemake@input[["own_mmms_signrand"]] 
))

print(list_own)

list_other <- unlist(list(
  snakemake@input[["other_sign_rand"]], 
  snakemake@input[["other_mark_rand"]], 
  snakemake@input[["other_mmms_rand"]], 
  snakemake@input[["other_mmms_signrand"]],
  snakemake@input[["other_mark_signrand"]] 
))

print(list_other)

#-------------------------------------------------------------------------------
# load data into lists and bind dfs together

df_own_list <- lapply(list_own, function(path){
  df <- base::readRDS(path)
  return(df)
})
df_own <- dplyr::bind_rows(df_own_list)

df_other_list <- lapply(list_other, function(path){
  df <- base::readRDS(path)
  return(df)
})
df_other <- dplyr::bind_rows(df_other_list)

#-------------------------------------------------------------------------------

# unify columns but later make it all in the original scripts

# for now, add an empty vector for resolution of own pval dfs
# but later, it should definitely be added

print(colnames(df_own))
print(colnames(df_other))

# match colnames to fit
df_own <- df_own[,match(colnames(df_other), colnames(df_own))]

#-------------------------------------------------------------------------------
# bind both dfs into one df, then correct pvals using BH = benjamini hochberg
# zebrafish p vals are excluded

df_all_pval <- dplyr::bind_rows(df_own, df_other)

df_all_pval$pval_corrected <- stats::p.adjust(df_all_pval$pval, method = "BH")

#-------------------------------------------------------------------------------

base::saveRDS(df_all_pval, snakemake@output[["pval_corrected_df"]])

utils::sessionInfo()
