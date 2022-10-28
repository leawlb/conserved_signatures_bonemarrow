#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(SingleR, quietly = TRUE)

sce <- readRDS(file = snakemake@input[["sce_06"]])

ref_all <- readRDS(file = snakemake@params[["ref_baccin_sce"]])
ref_hsc <- readRDS(file = snakemake@params[["ref_dahlin_sce"]])
ref_str <- readRDS(file = snakemake@params[["ref_dolgalev_sce"]])

# separate and fraction-specific annotation
sce_hsc <- sce[,sce$Fraction_ID == "hsc"]
sce_str <- sce[,sce$Fraction_ID == "str"]

pred_hsc <- SingleR(test=sce_hsc, ref=ref_hsc, 
                     labels=colData(ref_hsc)$identity_ref)
colData(sce_hsc)$Identity_ref_fraction <- pred_hsc$labels

pred_str <- SingleR(test=sce_str, ref=ref_str, 
                     labels=colData(ref_str)$identity_ref)
colData(sce_str)$Identity_ref_fraction <- pred_str$labels

# merge again and general annotation
sce <- cbind(sce_hsc, sce_str)

pred_all <- SingleR(test=sce, ref=ref_all, 
                    labels=colData(ref_all)$identity_ref)
colData(sce)$Identity_ref_all <- pred_all$labels
  
saveRDS(sce, file = snakemake@output[["sce_07"]])
