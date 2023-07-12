
library(SingleCellExperiment)

set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
sce_hsc <- readRDS(file = snakemake@input[["sce_hsc"]])
sce_str <- readRDS(file = snakemake@input[["sce_str"]])

# adjust SCE objects
# remove cells that are in the age merges but not the fraction merges
sce <- sce[,colnames(sce) %in% c(colnames(sce_hsc), colnames(sce_str))]

# add cell type information into "Identity" col
sce$Identity <- vector(length = ncol(sce))

# separate for faster calculation
sce_sep_str <- sce[,sce$Fraction_ID == "str"]
sce_sep_hsc <- sce[,sce$Fraction_ID == "hsc"]

# control that matches between different merges are perfect
stopifnot(!is.na(match(colnames(sce_sep_str), 
                       colnames(sce_str))))
table(is.na(match(colnames(sce_sep_str),
                  colnames(sce_str))))
table(is.na(match(colnames(sce_str), colnames(sce_sep_str))))

stopifnot(!is.na(match(colnames(sce_sep_hsc), 
                       colnames(sce_hsc))))
table(is.na(match(colnames(sce_sep_hsc), 
                  colnames(sce_str))))
table(is.na(match(colnames(sce_hsc), colnames(sce_sep_hsc))))

#-------------------------------------------------------------------------------

# add information: ADD CELL TYPE INFO INTO IDENTITY SLOT
# TODO: change to subcluster cell type once annotated
sce_sep_str$Identity <- sce_str$annotation_cluster[
  match(colnames(sce_sep_str), colnames(sce_str))]
sce_sep_hsc$Identity <- sce_hsc$annotation_cluster[
  match(colnames(sce_sep_hsc), colnames(sce_hsc))]

sce_output <- cbind(sce_sep_str, sce_sep_hsc)

print(colData(sce_output))

#-------------------------------------------------------------------------------

# add information: COPY IDs INTO ENSMUS_ID slot
rowData(sce_output)$ENSMUS_ID <- rowData(sce_output)$ID

#-------------------------------------------------------------------------------
saveRDS(sce_output, file = snakemake@output[["sce_output"]])