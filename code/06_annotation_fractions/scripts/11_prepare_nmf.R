#-------------------------------------------------------------------------------

library(DropletUtils)
library(tidyverse)
library(mclust)

#-------------------------------------------------------------------------------

# read corresponding sce object to trick snakemake into using the right wildcards
sce <- readRDS(file = snakemake@input[["sce_10"]])

fraction_curr <- sce$Fraction_ID[1]

factor_identities <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/factorCorrespondanceKitFrac.csv")
colnames(factor_identities)[1] <- "B6"

# NMF files provided by Adrien
# deposited there but should not change
if(fraction_curr == "hsc"){
  mcar_score <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/CarKit_clean.gene_spectra_score.k_26.dt_0_1.txt", sep = '\t')
  mcas_score <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/CastKit_clean.gene_spectra_score.k_20.dt_0_1.txt", sep = '\t')
  mmus_score <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/B6Kit_clean.gene_spectra_score.k_20.dt_0_1.txt", sep = '\t')
  mspr_score <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/SpretKit_clean.gene_spectra_score.k_20.dt_0_1.txt", sep = '\t')
  score_list <- list(mcar_score, mcas_score, mmus_score, mspr_score)
  names(score_list) <- c("Caro", "Cast", "B6", "Spret")
  species <- c("mcar", "mcas", "mmus", "mspr")
  
  mcar_usage <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/CarKit_clean.usages.k_26.dt_0_1.consensus.txt", sep = '\t')
  mcas_usage <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/CastKit_clean.usages.k_20.dt_0_1.consensus.txt", sep = '\t')
  mmus_usage <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/B6Kit_clean.usages.k_20.dt_0_1.consensus.txt", sep = '\t')
  mspr_usage <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/SpretKit_clean.usages.k_20.dt_0_1.consensus.txt", sep = '\t')
  usage_list<- list(mcar_usage, mcas_usage, mmus_usage, mspr_usage)
  names(usage_list) <- c("Caro", "Cast", "B6", "Spret")
}else if(fraction_curr == "str"){
  mmus_score <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/B6Str_clean.gene_spectra_score.k_19.dt_0_1.txt", sep = '\t')
  score_list <- list(mmus_score)
  names(score_list) <- c("B6")
  
  mmus_usage <- read.csv("/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/11_annotation/nmf_files/B6Str_clean.usages.k_19.dt_0_1.consensus.txt", sep = '\t')
  usage_list <- list(mmus_usage)
  names(usage_list) <- c("B6")
}

# for loop since multiple vectors/lists/names are used
score_list_temp <- score_list
usage_list_temp <- usage_list

for(i in 1:length(score_list_temp)){
  print(i)
  factor_identities_temp <- factor_identities
  factor_identities_temp <- factor_identities_temp[,c(2, which(colnames(factor_identities) == names(score_list_temp)[i]))]
  factor_identities_temp <- factor_identities_temp[!is.na(factor_identities_temp[,2]),]
  
  programs <- unlist(factor_identities_temp[2])
  print(programs)
  
  # change the score rownames to incorporate the 
  rownames(score_list_temp[[i]])[programs] <- factor_identities_temp$B6_factor_identity
  rownames(score_list_temp[[i]])[-programs] <- paste(fraction_curr,
                                                "program",
                                                rownames(score_list_temp[[i]])[-programs],
                                                species[i], sep = "_")
  rownames(score_list_temp[[i]]) <- sub("B6", "mmus", rownames(score_list_temp[[i]]))
  rownames(score_list_temp[[i]]) <- sub(" ", "_", rownames(score_list_temp[[i]]))
  print(rownames(score_list_temp[[i]]))
  
  colnames(score_list_temp[[i]])[colnames(score_list_temp[[i]]) == "X"] <- "ORIGINAL_PROGRAM"
  score_list_temp[[i]] <- as.data.frame(t(as.matrix(score_list_temp[[i]])))
  print(head(score_list_temp[[i]]))
  
  # change the usage colnames
  usage_list_temp[[i]] <- column_to_rownames(usage_list_temp[[i]], var = "X")
  colnames(usage_list_temp[[i]]) <- gsub("X", "", colnames(usage_list_temp[[i]]))
  colnames(usage_list_temp[[i]])[programs] <- factor_identities_temp$B6_factor_identity
  colnames(usage_list_temp[[i]])[-programs] <- paste(fraction_curr, 
                                                     "program",
                                                     colnames(usage_list_temp[[i]])[-programs], 
                                                     species[i], sep = "_")
  colnames(usage_list_temp[[i]]) <- sub("B6", "mmus", colnames(usage_list_temp[[i]]))
  colnames(usage_list_temp[[i]]) <- sub(" ", "_", colnames(usage_list_temp[[i]]))
  print(colnames(usage_list_temp[[i]]))
}

metadata(sce)$score_list <- score_list_temp
metadata(sce)$usage_list <- usage_list_temp

#-------------------------------------------------------------------------------

# cluster cells based on program expression, for each program

metadata(sce)$cluster_list <- metadata(sce)$usage_list

for(i in 1:length(metadata(sce)$cluster_list)){
  for(j in 1:ncol(metadata(sce)$cluster_list[[i]])){
    result <- Mclust(metadata(sce)$cluster_list[[i]][,j], G = 2)
    metadata(sce)$cluster_list[[i]][,j] <- result$classification
  }
}

names(metadata(sce)$score_list) <- species
names(metadata(sce)$usage_list) <- species
names(metadata(sce)$cluster_list) <- species

saveRDS(sce, file = snakemake@output[["sce_11"]])

