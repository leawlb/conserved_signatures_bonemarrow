
library(SingleCellExperiment)
set.seed(37)

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_input"]])
min_cells <- snakemake@params[["min_cells"]]

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

cell_nrs <- table(sce$Identity)
print(min_cells)
for(i in 1:length(cell_nrs)){
  if(cell_nrs[i] <= min_cells){
    print(paste("removing", names(cell_nrs)[i]))
    sce <- sce[,!sce$Identity == names(cell_nrs)[i]]
  }
}

print(levels(sce$Identity))
sce$Identity <- unfactor(sce$Identity)
print(assignment_list$Identity[assignment_list$Identity %in% sce$Identity])

sce$Identity <- factor(sce$Identity, levels = assignment_list$Identity[
  assignment_list$Identity %in% sce$Identity])
 print(levels(sce$Identity))
 
print(colData(sce))

#-------------------------------------------------------------------------------
saveRDS(sce, file = snakemake@output[["sce_output"]])