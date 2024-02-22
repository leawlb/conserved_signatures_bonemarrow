library(SingleCellExperiment)
set.seed(37)

source(snakemake@params[["metrics_functions"]])

sce <- readRDS(file = snakemake@input[["sce_input"]])
cci <- readRDS(file = snakemake@input[["cci_input"]])

### STEP 1: IDENTIY PAIR METRICS
ipi_list <- extract_ident_pair_info(cci = cci, sce = sce)

### STEP 2: IDENTITY METRICS
idi_list <- extract_ident_info(ipi_list = ipi_list)

### STEP 3 AND 4: IDENTITY LIGAND/RECEPTOR INFO
ident_lrs_list <- lapply(idi_list[[2]], extract_lrs_info)

ident_nrlrs <- extract_lrs_nrs(cci = cci, lrs_list = ident_lrs_list)

ilr_list <- list("overview" = ident_nrlrs, "separate" = ident_lrs_list)

### save
metrics_lists <- list("ipi" = ipi_list, "idi" = idi_list, "ilr" = ilr_list)
saveRDS(metrics_lists, file = snakemake@output[["metrics_lists"]])