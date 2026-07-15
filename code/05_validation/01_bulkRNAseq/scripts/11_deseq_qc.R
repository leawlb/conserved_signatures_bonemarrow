
# We want to analyse within target populations across conditions, so separate
# into two objects, one for Neo1, and one for Trkb samples

# perform several analyses on dds objects for QC purposes
# then transform for DESeq2 analysis and export analysis

#-------------------------------------------------------------------------------

library(DESeq2)
library(tidyverse)

# load
input_path <- snakemake@input[["input_path_dds"]]

dds <- readRDS(input_path)

fc_cutoff <- 0.25 # 
alpha_cutoff <- 0.1 # default

ab_target <- snakemake@wildcards[["ab_target"]] 

#-------------------------------------------------------------------------------

# subset to required samples by target
dds_subset <- dds[,dds$Target_ID == ab_target]

# mmus_msc_trkb_p_08 has consistently been low quality
# both in multiQC (lowest numbers of reads) and in vsd results when included
# and with mmus_msc_trkb_p_08 included there were no significant results from DESeq 

# also change design; we want to compare between pos and neg conditions

print(colData(dds_subset))

if(ab_target == "trkb"){
  dds_subset <- dds_subset[,!dds_subset$Object_ID == "mmus_msc_trkb_p_08"]
  DESeq2::design(dds_subset) <- ~ Mouse_ID + Condition_ID

}else if(ab_target == "neo1"){
  DESeq2::design(dds_subset) <- ~ Mouse_ID + Condition_ID
}

dds_subset$Object_ID <- factor(
  dds_subset$Object_ID,
  levels = levels(dds_subset$Object_ID)[
  levels(dds_subset$Object_ID) %in% unique(unfactor(dds_subset$Object_ID))])

# filter out almost empty rows (non-informative genes) across 6 samples, each
dds_subset <- dds_subset[rowSums(counts(dds_subset)) > 20,]

# set size factors:
dds_subset <- estimateSizeFactors(dds_subset)
sizeFactors(dds_subset)

#-------------------------------------------------------------------------------
# transform for expression value visualisation only
# using blind = FALSE assuming the design is suitable after qc and adjustments

vsd <- DESeq2::vst(dds_subset, blind = FALSE)

vsd_output_path <- snakemake@output[["output_path_vsd"]]
saveRDS(vsd, vsd_output_path)

#-------------------------------------------------------------------------------
#transform and get results 

tdsq <- DESeq2::DESeq(dds_subset)

res <- DESeq2::results(
  tdsq,
  lfcThreshold=fc_cutoff, 
  contrast=c("Condition_ID", "p", "n"),
  alpha = alpha_cutoff)

results_df <- res %>%
  as.data.frame() %>%
  arrange(desc(abs(log2FoldChange)))

results_df[1:100,]

results_output_path <- snakemake@output[["output_path_results"]]
saveRDS(results_df, results_output_path)

sessionInfo()
