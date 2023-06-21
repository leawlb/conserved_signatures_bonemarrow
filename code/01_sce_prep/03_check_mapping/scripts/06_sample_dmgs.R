#-------------------------------------------------------------------------------

library(SingleCellExperiment, quietly = TRUE) 
library(Seurat, quietly = TRUE)
library(biomaRt, quietly = TRUE)
library(scran, quietly = TRUE)

# load og (one genome) and fg (four genomes) objects
sce_og <- readRDS(file = snakemake@input[["sce_02"]])
sce_fg <- readRDS(file = snakemake@input[["sce_fg"]])

nr_hvgs <- snakemake@params[["nr_hvgs"]]
logFC_sample_mgs <- snakemake@params[["logFC_sample_mgs"]]
minPC_sample_mgs <- snakemake@params[["minPC_sample_mgs"]]
  
sce_fg <- sce_fg[,which(!is.na(match(colnames(sce_fg), colnames(sce_og))))]

print(sce_og)
print(sce_fg)

#-------------------------------------------------------------------------------

name_curr <- colData(sce_og)$Object_ID[1] # name to display in plots
# add symbols to mspr objects 
# This part was taken from Perrine's script on how to download the relevant data 
if(grepl("mspr", name_curr)){

  ensembl_mspr <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                          host="https://www.ensembl.org",
                          dataset="mspretus_gene_ensembl")
  ensembl_list_mspr=getBM(attributes=c("ensembl_gene_id",
                                       "mmusculus_homolog_ensembl_gene",
                                       "mmusculus_homolog_associated_gene_name"), 
                          mart=ensembl_mspr)

  intersect_IDs <- intersect(rowData(sce_fg)$Symbol, 
                             ensembl_list_mspr$ensembl_gene_id)

  ensembl_list_mspr <- ensembl_list_mspr[
    ensembl_list_mspr$ensembl_gene_id %in% intersect_IDs,]
  sce_fg <- sce_fg[rowData(sce_fg)$Symbol %in% intersect_IDs,]

  ensembl_list_mspr <- ensembl_list_mspr[
    !duplicated(ensembl_list_mspr$ensembl_gene_id),]
  sce_fg <- sce_fg[!duplicated(rowData(sce_fg)$Symbol),]

  rowData(sce_fg)$Symbol[
    match(
      ensembl_list_mspr$ensembl_gene_id,
      rownames(sce_fg))] <- ensembl_list_mspr$mmusculus_homolog_associated_gene_name

  rowData(sce_fg)$Symbol[
    rowData(sce_fg)$Symbol == ""] <- rowData(sce_fg)$ID[
      rowData(sce_fg)$Symbol == ""]
}

#-------------------------------------------------------------------------------

# make og and fg ready for comparison

# rows
shared_genes <- intersect(rowData(sce_og)$Symbol, rowData(sce_fg)$Symbol)

sce_og_shared <- sce_og[rowData(sce_og)$Symbol %in% shared_genes,]
sce_fg_shared <- sce_fg[rowData(sce_fg)$Symbol %in% shared_genes,]
sce_og_shared <- sce_og_shared[order(rowData(sce_og_shared)$Symbol),]
sce_fg_shared <- sce_fg_shared[order(rowData(sce_fg_shared)$Symbol),]
sce_og_shared <- sce_og_shared[!duplicated(rowData(sce_og_shared)$Symbol),]
sce_fg_shared <- sce_fg_shared[!duplicated(rowData(sce_fg_shared)$Symbol),]

# columns
colnames(sce_og_shared) <- paste0(colnames(sce_og_shared), "_og")
sce_og_shared$Genome <- rep("OneGenome", ncol(sce_og_shared))

colnames(sce_fg_shared) <- paste0(colnames(sce_fg_shared), "_fg")
sce_fg_shared$Genome <- rep("FourGenomes", ncol(sce_fg_shared))

colnames(rowData(sce_og_shared))[colnames(rowData(sce_og_shared)) == "ID"] <- "ID_og"
colnames(rowData(sce_fg_shared))[colnames(rowData(sce_fg_shared)) == "ID"] <- "ID_fg"

rownames(sce_og_shared) <- rowData(sce_og_shared)$Symbol 
rownames(sce_fg_shared) <- rowData(sce_fg_shared)$Symbol 

colData(sce_og_shared) <- colData(sce_og_shared)[,-20]

# combine
sce_tog <- cbind(sce_og_shared, sce_fg_shared)
sce_tog

#-------------------------------------------------------------------------------

quick_clust <- quickCluster(sce_tog)
sce_tog <- computeSumFactors(sce_tog, cluster = quick_clust)
sce_tog <- logNormCounts(sce_tog) 

genevar <- modelGeneVar(sce_tog)
hvg <- getTopHVGs(genevar, n=nr_hvgs)

#-------------------------------------------------------------------------------

# get markergenes = genes that are differentially "expressed" (=mapped)

# convert to seurat
seurat <- as.Seurat(sce_tog, counts = "counts", data = "logcounts",) 
Idents(seurat) <- sce_tog$Genome

# get hvgs for marker genes
hvgs <- modelGeneVar(sce_tog)
hvgs <- getTopHVGs(hvgs, n=nr_hvgs)

# find marker genes for each genome
genomes <- unique(sce_tog$Genome)
clustlist <- as.list(genomes)
cluster_markers <- lapply(clustlist, function(x){
  markers <- FindMarkers(seurat, test.use = "wilcox", ident.1 = x, 
                         features = hvgs, logfc.threshold = logFC_sample_mgs,
                         min.pct = minPC_sample_mgs)
  markers <- markers[order(abs(markers$avg_log2FC), decreasing=TRUE),]
  markers$which_cluster <- rep(x, nrow(markers))
  return(markers)
})

names(cluster_markers) <- genomes
markergenes_og <- cluster_markers[["OneGenome"]]

saveRDS(markergenes_og, snakemake@output[["dmgs"]])