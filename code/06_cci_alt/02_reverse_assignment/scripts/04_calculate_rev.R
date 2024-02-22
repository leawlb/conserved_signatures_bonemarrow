library(SingleCellExperiment)
set.seed(37)

source(snakemake@params[["main_functions"]])
source(snakemake@params[["alternative_functions"]])

sce <- readRDS(file = snakemake@input[["sce_input"]])
lrdb <- readRDS(file = snakemake@input[["lrdb"]])

min_perc <- snakemake@params[["min_perc"]]
top_level <- snakemake@params[["top_level"]]
assay_use <- snakemake@params[["assay_use"]]
condition <- snakemake@params[["condition"]]

print(paste("cutoff is", min_perc))
print(paste("top_level is", top_level))
print(paste("assay_use is", assay_use))

if(top_level == 0){top_level <- NULL}

cci <- calculate_scores_wrap(sce = sce, lrdb = lrdb, 
                             cutoff = min_perc, 
                             top_level = top_level, 
                             assay_use = assay_use,
                             condition = condition)

cci$Identities$Species <- rep(sce$Species_ID[1], nrow(cci$Identities))
cci$Identities$Age <- rep(sce$Age_ID[1], nrow(cci$Identities))

saveRDS(cci, snakemake@output[["cci_output"]])