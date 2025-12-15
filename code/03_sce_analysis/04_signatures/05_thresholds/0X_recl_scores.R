
# get reclustering scores for reclustered datasets from other species

# determine random number generator for sample
library(parallel)
RNGkind("L'Ecuyer-CMRG") # using this for reclustering

#-------------------------------------------------------------------------------

library(mclust, quietly = TRUE)
library(bluster, quietly = TRUE)
library(S4Vectors, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(SeuratObject, quietly = TRUE)

source(file = snakemake@params[["reclustering_functions"]])

#-------------------------------------------------------------------------------
# load seurat object lists with reclustered labels
# ordered by the different thresholds on analyses used to generate 
# signature gene lists

seu_list_all <- base::readRDS(snakemake@input[["seu_list"]])
#seu_list_all <- readRDS("/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/04_signatures/05_thresholds/ts_hscs_progenitors_recl_by_p.rds")

threshold_used <- snakemake@params[["threshold_used"]]

#-------------------------------------------------------------------------------
# for each item of the list, add the used to misc so it's recorded for later

for(threshold in names(seu_list_all)){
  seu_list_all[[threshold]]@misc$threshold_used <- threshold
}

if(threshold_used == "conserved_signature_treshold_testing_by_t"){
  # below a certain number of genes there are not enough pcs to calculate scores
  seu_list_all <- seu_list_all[!names(seu_list_all) %in% c("0.7", "1")]
  names(seu_list_all)
}

#-------------------------------------------------------------------------------
# for each item of the list, calculate the reclustering scores 
score_df_list <- lapply(seu_list_all, function(seu){
    
    print(seu@misc$used_genes)
    print(seu@misc$threshold_used)

    score_df <- calculate_scores(seu) # own function
    
    # add specific info
    score_df$conservation_level <- base::rep(seu@misc$used_genes,
                                             nrow(score_df))
    score_df$resolution <- base::rep(seu@misc$resolution,
                                             nrow(score_df))   
    score_df$nr_genes_used <- base::rep(seu@misc$nr_genes_used,
                                        nrow(score_df))
    score_df$threshold_used <- base::rep(seu@misc$threshold_used,
                                        nrow(score_df))
    
    return(score_df)
})

names(score_df_list) <- names(seu_list_all)

#-------------------------------------------------------------------------------

base::saveRDS(score_df_list, snakemake@output[["score_df_list"]])

utils::sessionInfo()
