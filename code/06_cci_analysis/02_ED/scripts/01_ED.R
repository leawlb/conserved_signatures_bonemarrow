set.seed(37)

cci <- readRDS(file = snakemake@input[["cci_input"]]) # cci prep output

#scaling 
cci_scaled <- scale(t(cci$Score))

#euclidian distances
cci_dist <- dist(cci_scaled, method = 'euclidean')

# make a df for convenience
cci_mat <- as.matrix(cci_dist)
cci_dist_df <- as.data.frame(cci_mat)

saveRDS(cci_dist_df, file = snakemake@output[["cci_dist_output"]])
