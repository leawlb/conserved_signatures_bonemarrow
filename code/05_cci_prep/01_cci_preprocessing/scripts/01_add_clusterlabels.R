
library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
sce_hsc <- readRDS(file = snakemake@input[["sce_hsc"]])
sce_str <- readRDS(file = snakemake@input[["sce_str"]])

coln_ident <- snakemake@params[["colname_identity"]]
print(coln_ident)
# making this modifiable in config causes sce$col to turn into colData(sce)[,colnames(colData(sce)) == col]
print(unique(colData(sce_hsc)[,colnames(colData(sce_hsc)) == coln_ident]))
print(unique(colData(sce_str)[,colnames(colData(sce_str)) == coln_ident]))

stopifnot(!is.na(unique(colData(sce_hsc)[,colnames(colData(sce_hsc)) == coln_ident])))
stopifnot(!is.na(unique(colData(sce_str)[,colnames(colData(sce_str)) == coln_ident])))

#-------------------------------------------------------------------------------

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

# add information from stated coldata slot into IDENTITY slot
print(levels(colData(sce_str)[,colnames(colData(sce_str)) == coln_ident]))
print(levels(colData(sce_hsc)[,colnames(colData(sce_hsc)) == coln_ident]))

sce_sep_str$Identity <- colData(sce_str)[,colnames(colData(sce_str)) == coln_ident][
  match(colnames(sce_sep_str), colnames(sce_str))]
sce_sep_hsc$Identity <- colData(sce_hsc)[,colnames(colData(sce_hsc)) == coln_ident][
  match(colnames(sce_sep_hsc), colnames(sce_hsc))]

sce_output <- cbind(sce_sep_str, sce_sep_hsc)

print(colData(sce_output))

#-------------------------------------------------------------------------------

# add information: COPY IDs INTO ENSMUS_ID slot
rowData(sce_output)$ENSMUS_ID <- rowData(sce_output)$ID

#-------------------------------------------------------------------------------
saveRDS(sce_output, file = snakemake@output[["sce_output"]])