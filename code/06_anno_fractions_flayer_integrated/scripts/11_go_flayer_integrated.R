#-------------------------------------------------------------------------------

library(DropletUtils)
library(scran)
library(GOfuncR)
set.seed(37)

cluster_markers <- readRDS(file = snakemake@input[["results_markers"]] )

# use marker genes as input for GO analysis
get_go_results <- function(markers){
  
  go_input_df <- data.frame(
    gene_ids = rownames(markers),
    gene_scores =  markers$avg_log2FC
  )
  
  # perform GO analysis using logFC to determine enrichment
  go_result <- go_enrich(go_input_df, organismDb='Mus.musculus',
                         test='wilcoxon', n_randsets=100)
  return(go_result)  
}

go_results <- lapply(cluster_markers, get_go_results)

saveRDS(go_results, file = snakemake@output[["results_go"]] )

