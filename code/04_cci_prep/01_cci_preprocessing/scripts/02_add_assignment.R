
library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
min_cells <- snakemake@params[["min_cells"]]
min_cells <- 80

assignment_list <- read.csv(snakemake@input[["assignment"]], 
                       header = TRUE, 
                       sep = ";", 
                       check.names=FALSE, 
                       stringsAsFactors=FALSE, 
                       as.is=TRUE, 
                       colClasses = "character")

print(assignment_list)

stopifnot(unique(sce$Identity) %in% unique(assignment_list$Identity))
stopifnot(unique(assignment_list$Identity) %in% unique(sce$Identity))

sce$Assignment <- vector(length = ncol(sce))
for(i in unique(sce$Identity)){
  sce$Assignment[sce$Identity == i] <- assignment_list$Assignment[
    assignment_list$Identity == i]
}

print(ncol(sce))
sce <- sce[,!sce$Assignment == "remove"]
print(ncol(sce))

# remove identities with low cell numbers
sce$Identity <- unfactor(sce$Identity)
ident_rm <- names(table(sce$Identity))[which(table(sce$Identity) <= min_cells)]

print(paste("removing", ident_rm))
# only if the ident_rm to be removed are still in SCE object
sce <- sce[,!sce$Identity %in% ident_rm]

# factor again using assignment_list for the correct sequence of identities
print(assignment_list$Identity[assignment_list$Identity %in% sce$Identity])

sce$Identity <- factor(sce$Identity, levels = assignment_list$Identity[
  assignment_list$Identity %in% sce$Identity])
 print(levels(sce$Identity))
 
print(colData(sce))

#-------------------------------------------------------------------------------
saveRDS(sce, file = snakemake@output[["sce_output"]])