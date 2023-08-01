set.seed(37)
library(stringr)

#-------------------------------------------------------------------------------

cci <- readRDS(file = snakemake@input[["cci_input"]]) # cci prep output

# prepare data for PCA calculation
temp_df <- data.frame(t(cci$Score))
temp_df[is.na(temp_df)] <- 0
temp_df <- temp_df[,-which(colSums(temp_df) == 0)]

print(temp_df[1:10,1:10])

# use base::pcromp for PCA calculation
pca_df <- prcomp(temp_df, scale. = TRUE)

str(pca_df)

# add info to df more suited for visualisation
ggdf <- as.data.frame(pca_df$x[,1:5])

ggdf$Age <- vector(length = nrow(ggdf))
ggdf$Species <- vector(length = nrow(ggdf))
ggdf$ident_pair <- vector(length = nrow(ggdf))
ggdf$emitter <- vector(length = nrow(ggdf))
ggdf$receiver <- vector(length = nrow(ggdf))

for(i in 1:nrow(ggdf)){
  ggdf$ident_pair[i] <- as.character(str_split(rownames(ggdf)[i], "_")[[1]][1])
  ggdf$Age[i] <- as.character(str_split(rownames(ggdf)[i], "_")[[1]][3])
  ggdf$Species[i] <- as.character(str_split(rownames(ggdf)[i], "_")[[1]][2])
  
  ggdf$emitter[i] <- as.character(str_split(ggdf$ident_pair[i], "&")[[1]][1])
  ggdf$receiver[i] <- as.character(str_split(ggdf$ident_pair[i], "&")[[1]][2])
  
}

# export both
saveRDS(list(pca_df, ggdf), file = snakemake@output[["pca_output"]])