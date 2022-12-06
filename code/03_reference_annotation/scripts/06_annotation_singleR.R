#-------------------------------------------------------------------------------

library(DropletUtils, quietly = TRUE)
library(SingleR, quietly = TRUE)

sce <- readRDS(file = snakemake@input[["sce_04"]])

ref_baccin <- readRDS(file = snakemake@params[["ref_baccin_sce"]])
ref_dahlin <- readRDS(file = snakemake@params[["ref_dahlin_sce"]])
ref_dolgalev <- readRDS(file = snakemake@params[["ref_dolgalev_sce"]])
ref_lipka <- readRDS(file = snakemake@params[["ref_lipka_sce"]])

# separate and fraction-specific annotation
sce_hsc <- sce[,sce$Fraction_ID == "hsc"]
sce_str <- sce[,sce$Fraction_ID == "str"]

pred_hsc1 <- SingleR(test=sce_hsc, ref=ref_dahlin, 
                    labels=colData(ref_dahlin)$identity_ref)
colData(sce_hsc)$Identity_ref_fraction1 <- pred_hsc1$labels
pred_hsc2 <- SingleR(test=sce_hsc, ref=ref_lipka, 
                     labels=colData(ref_lipka)$identity_ref)
colData(sce_hsc)$Identity_ref_fraction2 <- pred_hsc2$labels
print("... finished hscs")

pred_str1 <- SingleR(test=sce_str, ref=ref_dolgalev, 
                    labels=colData(ref_dolgalev)$identity_ref)
colData(sce_str)$Identity_ref_fraction1 <- pred_str1$labels
pred_str2 <- SingleR(test=sce_str, ref=ref_dolgalev, 
                     labels=colData(ref_dolgalev)$labelsimple)
colData(sce_str)$Identity_ref_fraction2 <- pred_str2$labels
print("... finished stromal")

# merge again and general annotation
sce <- cbind(sce_hsc, sce_str)

pred_all <- SingleR(test=sce, ref=ref_baccin, 
                    labels=colData(ref_baccin)$identity_ref)
colData(sce)$Identity_ref_all <- pred_all$labels
print(colData(sce))

saveRDS(sce, file = snakemake@output[["sce_06"]])
