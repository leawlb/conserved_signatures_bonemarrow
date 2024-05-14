#-------------------------------------------------------------------------------
# basic GO analysis

#library(scran)
library(GOfuncR)
set.seed(37)

#-------------------------------------------------------------------------------

cluster_markers <- readRDS(file = snakemake@input[["markers"]])

#-------------------------------------------------------------------------------
# use marker genes as input for GO analysis
get_go_results <- function(markers){
  
  go_input_df <- data.frame(
    gene_ids = rownames(markers),
    gene_scores =  markers$avg_log2FC
  )
  if(nrow(go_input_df) > 50){
    go_input_df <- go_input_df[1:50,] 
  }

  # perform GO analysis using logFC to determine enrichment
  go_result <- GOfuncR::go_enrich(go_input_df, organismDb='Mus.musculus',
                                  test='wilcoxon', n_randsets=100)
  return(go_result)  
}

go_results <- lapply(cluster_markers, get_go_results)

#-------------------------------------------------------------------------------
saveRDS(go_results, file = snakemake@output[["go"]])

sessionInfo()
