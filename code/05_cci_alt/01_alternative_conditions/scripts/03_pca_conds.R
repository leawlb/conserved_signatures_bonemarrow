set.seed(37)
library(stringr)

#-------------------------------------------------------------------------------

cci <- readRDS(file = snakemake@input[["cci_input"]]) # cci prep output

# prepare data for PCA calculation
temp_df <- data.frame(t(cci$Score))
temp_df[is.na(temp_df)] <- 0
temp_df <- temp_df[,-which(colSums(temp_df) == 0)]

print(temp_df[1:10,1:10])

# use base::prcomp for PCA calculation
pca_df <- prcomp(temp_df, scale. = TRUE)

str(pca_df)

# add info to a df more suited for visualisation
ggdf <- as.data.frame(pca_df$x[,1:5])

ggdf$ident_pair <- rownames(ggdf)
ggdf$age <- rep(cci$Identities$age[1], nrow(ggdf))
ggdf$species <- rep(cci$Identities$species[1], nrow(ggdf))
ggdf$condition <- rep(cci$Identities$condition[1], nrow(ggdf))

ggdf$emitter <- vector(length = nrow(ggdf))
ggdf$receiver <- vector(length = nrow(ggdf))

for(i in 1:nrow(ggdf)){
  ggdf$emitter[i] <- as.character(str_split(ggdf$ident_pair[i], "&")[[1]][1])
  ggdf$receiver[i] <- as.character(str_split(ggdf$ident_pair[i], "&")[[1]][2])
  
}

sum_df <- summary(pca_df)

PC1_var <- sum_df$importance[2,1]*100
PC2_var <- sum_df$importance[2,2]*100
PC3_var <- sum_df$importance[2,3]*100
PC4_var <- sum_df$importance[2,4]*100
PC5_var <- sum_df$importance[2,5]*100

ggdf$PC1_var <- PC1_var
ggdf$PC2_var <- PC2_var
ggdf$PC3_var <- PC3_var
ggdf$PC4_var <- PC4_var
ggdf$PC5_var <- PC5_var

print(colnames(ggdf))

saveRDS(list(pca_df, ggdf), file = snakemake@output[["pca_output"]])