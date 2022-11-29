# Prepare NMF objects and add to SCE for further use
#-------------------------------------------------------------------------------
library(DropletUtils)
library(tidyverse)
#-------------------------------------------------------------------------------
sce_07 <- readRDS(snakemake@input[["sce_07"]]) 
nmf_path <- snakemake@input[["nmf"]] 
nr_hvgs <- snakemake@params[["nr_hvgs"]]
source(file = snakemake@params[["sce_functions"]])
#-------------------------------------------------------------------------------

# make sure that the correct object is loaded into the correct variable 
print("loading hsc_usages")
print(nmf_path[1])
hsc_usages <- read.csv(nmf_path[1], sep = '\t')

print("loading hsc_scores")
print(nmf_path[2])
hsc_scores <- read.csv(nmf_path[2], sep = '\t')

print("loading str_usages")
print(nmf_path[3])
str_usages <- read.csv(nmf_path[3], sep = '\t')

print("loading str_scores")
print(nmf_path[4])
str_scores <- read.csv(nmf_path[4], sep = '\t')

#-------------------------------------------------------------------------------
# prepare all nmf objects
# HSC
hsc_usages <- column_to_rownames(hsc_usages, var = "X")
colnames(hsc_usages) <- gsub("X", "hsc_program_", colnames(hsc_usages))
hsc_usages <- hsc_usages + 1
hsc_usages <- log10(hsc_usages)

hsc_scores <- column_to_rownames(hsc_scores, var = "X")
rownames(hsc_scores) <- paste0("hsc_program_", rownames(hsc_scores))
hsc_scores <- as.data.frame(t(as.matrix(hsc_scores)))

# STR
str_usages <- column_to_rownames(str_usages, var = "X")
colnames(str_usages) <- gsub("X", "str_program_", colnames(str_usages))
str_usages <- str_usages + 1
str_usages <- log10(str_usages)

str_scores <- column_to_rownames(str_scores, var = "X")
rownames(str_scores) <- paste0("str_program_", rownames(str_scores))
str_scores <- as.data.frame(t(as.matrix(str_scores)))

#-------------------------------------------------------------------------------
# add gene scores as metadata

metadata(sce_07)$Genescores_NMF_hsc <- hsc_scores
metadata(sce_07)$Genescores_NMF_str <- str_scores

#-------------------------------------------------------------------------------
# add barcode program usages as colData
# this requires separation of stromal and hsc programs

# add hsc programs to hsc positions as well as empty slots in str positions
hsc_prog_hsc <- hsc_usages
hsc_prog_str <- data.frame(row.names = rownames(colData(sce_07)[is.na(match(
  rownames(colData(sce_07)), rownames(hsc_prog_hsc))),])) 
print(hsc_prog_str)
hsc_prog_str[,1:ncol(hsc_prog_hsc)] <- rep(0, nrow(hsc_prog_str))
colnames(hsc_prog_str) <- colnames(hsc_prog_hsc)
hsc_prog <- rbind(hsc_prog_hsc, hsc_prog_str)

colData(sce_07) <- cbind(colData(sce_07), hsc_prog[match(
  rownames(colData(sce_07)), rownames(hsc_prog)),])

# add stromal programs to str positions as well as empty slots in hsc positions
str_prog_str <- str_usages
str_prog_hsc <- data.frame(row.names = rownames(colData(sce_07)[is.na(match(
  rownames(colData(sce_07)), rownames(str_prog_str))),])) 
str_prog_hsc[,1:ncol(str_prog_str)] <- rep(0, nrow(str_prog_hsc))
colnames(str_prog_hsc) <- colnames(str_prog_str)
str_prog <- rbind(str_prog_str, str_prog_hsc)

colData(sce_07) <- cbind(colData(sce_07), str_prog[match(
  rownames(colData(sce_07)), rownames(str_prog)),])

#-------------------------------------------------------------------------------
# dimensionality reduction to add to the saved objects
sce_07 <- reduce_dims(sce_07, nr_hvgs = nr_hvgs) # own function
print(sce_07)

#-------------------------------------------------------------------------------
# save object containing new info
saveRDS(sce_07, snakemake@output[["sce_08"]]) 