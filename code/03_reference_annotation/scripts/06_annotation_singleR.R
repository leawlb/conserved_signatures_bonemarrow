#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(SingleR, quietly = TRUE)
set.seed(37)

sce <- readRDS(file = snakemake@input[["sce_05"]])

ref_baccin <- readRDS(file = snakemake@params[["ref_baccin_sce"]])
ref_dahlin <- readRDS(file = snakemake@params[["ref_dahlin_sce"]])
ref_dolgalev <- readRDS(file = snakemake@params[["ref_dolgalev_sce"]])

# separate and fraction-specific annotation
sce_hsc <- sce[,sce$Fraction_ID == "hsc"]
sce_str <- sce[,sce$Fraction_ID == "str"]

pred_hsc <- SingleR(test=sce_hsc, ref=ref_dahlin, 
                    labels=colData(ref_dahlin)$identity_ref)
colData(sce_hsc)$Identity_ref_fraction <- pred_hsc$labels
saveRDS(pred_hsc, file = snakemake@output[["pred_hsc"]])
print("... finished hscs")

pred_str <- SingleR(test=sce_str, ref=ref_dolgalev, 
                    labels=colData(ref_dolgalev)$labelsimple)
colData(sce_str)$Identity_ref_fraction <- pred_str$labels
saveRDS(pred_str, file = snakemake@output[["pred_str"]])
print("... finished stromal")

# merge again and general annotation
sce <- cbind(sce_hsc, sce_str)

pred_all <- SingleR(test=sce, ref=ref_baccin, 
                    labels=colData(ref_baccin)$identity_ref)
colData(sce)$Identity_ref_all <- pred_all$labels
saveRDS(pred_all, file = snakemake@output[["pred_all"]])
print(colData(sce))

saveRDS(sce, file = snakemake@output[["sce_06"]])
