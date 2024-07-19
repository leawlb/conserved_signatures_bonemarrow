#-------------------------------------------------------------------------------
# basic GO analysis only

library(GOfuncR)
set.seed(37)

#-------------------------------------------------------------------------------

cluster_markers <- base::readRDS(file = snakemake@input[["markers"]])

#-------------------------------------------------------------------------------
# use marker genes as input for GO analysis
get_go_results <- function(markers){
  
  go_input_df <- base::data.frame(
    gene_ids = rownames(markers),
    gene_scores =  markers$avg_log2FC
  )
  if(nrow(go_input_df) > 50){
    go_input_df <- go_input_df[1:50,] 
  }

  # perform GO analysis with wilcoxon test on avg_log2FC to determine enrichment
  go_result <- GOfuncR::go_enrich(
    go_input_df, 
    organismDb='Mus.musculus', # using MM since it's only a basic analysis 
    test='wilcoxon', 
    n_randsets=100
    )
  
  return(go_result)  
}

go_results <- lapply(cluster_markers, get_go_results)

#-------------------------------------------------------------------------------

base::saveRDS(go_results, file = snakemake@output[["go"]])

utils::sessionInfo()