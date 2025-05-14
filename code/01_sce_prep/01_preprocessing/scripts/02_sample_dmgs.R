#-------------------------------------------------------------------------------
# This is a long script for one step: Finding differentially mapped genes (dmgs)
# between individual SCE objects aligned with different reference genomes.
# Prepare for direct comparison between both SCEs, then find marker genes.

library(SingleCellExperiment, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(scran, quietly = TRUE)

#-------------------------------------------------------------------------------
# load og (one genome) and fg (four genomes) objects
sce_og <- base::readRDS(file = snakemake@input[["sce_input"]])
sce_fg <- base::readRDS(file = snakemake@input[["sce_fg"]])

nr_hvgs <- snakemake@params[["nr_hvgs"]]
logFC_sample_dmgs <- snakemake@params[["logFC_sample_dmgs"]]

ensembl_list_mspr <- base::readRDS(snakemake@input[["ensembl_list_mspr"]])
head(ensembl_list_mspr)

sce_fg <- sce_fg[,which(!is.na(base::match(colnames(sce_fg), 
                                           colnames(sce_og))))]

print(sce_og)
print(sce_fg)

name_curr <- colData(sce_og)$Object_ID[1] # name to display in plots (later)

#-------------------------------------------------------------------------------

# add symbols to mspr fg objects, which are missing annotations
# ensembl list for mspr was previously downloaded and pre-processed in 
# metadata/scRNAseq/01_sce_prep/ensembl_mspr/
if(grepl("mspr", name_curr)){
  
  intersect_IDs <- BiocGenerics::intersect(
    rowData(sce_fg)$Symbol, ensembl_list_mspr$ensembl_gene_id)

  ensembl_list_mspr <- ensembl_list_mspr[
    ensembl_list_mspr$ensembl_gene_id %in% intersect_IDs,]
  sce_fg <- sce_fg[rowData(sce_fg)$Symbol %in% intersect_IDs,]

  ensembl_list_mspr <- ensembl_list_mspr[
    !base::duplicated(ensembl_list_mspr$ensembl_gene_id),]
  sce_fg <- sce_fg[!base::duplicated(rowData(sce_fg)$Symbol),]

  rowData(sce_fg)$Symbol[
    match(
      ensembl_list_mspr$ensembl_gene_id,
      rownames(sce_fg))] <- ensembl_list_mspr$mmusculus_homolog_associated_gene_name

  rowData(sce_fg)$Symbol[
    rowData(sce_fg)$Symbol == ""] <- rowData(sce_fg)$ID[
      rowData(sce_fg)$Symbol == ""]
}

#-------------------------------------------------------------------------------
# get og and fg ready for comparison

# rows
shared_genes <- BiocGenerics::intersect(rowData(sce_og)$Symbol, 
                                        rowData(sce_fg)$Symbol)

sce_og_shared <- sce_og[rowData(sce_og)$Symbol %in% shared_genes,]
sce_fg_shared <- sce_fg[rowData(sce_fg)$Symbol %in% shared_genes,]

sce_og_shared <- sce_og_shared[base::order(rowData(sce_og_shared)$Symbol),]
sce_fg_shared <- sce_fg_shared[base::order(rowData(sce_fg_shared)$Symbol),]

sce_og_shared <- sce_og_shared[
  !base::duplicated(rowData(sce_og_shared)$Symbol),]
sce_fg_shared <- sce_fg_shared[
  !base::duplicated(rowData(sce_fg_shared)$Symbol),]

rownames(sce_og_shared) <- rowData(sce_og_shared)$Symbol 
rownames(sce_fg_shared) <- rowData(sce_fg_shared)$Symbol 

# columns
colnames(sce_og_shared) <- base::paste0(colnames(sce_og_shared), "_og")
colnames(sce_fg_shared) <- base::paste0(colnames(sce_fg_shared), "_fg")

sce_og_shared$Genome <- base::rep("OneGenome", ncol(sce_og_shared))
sce_fg_shared$Genome <- base::rep("FourGenomes", ncol(sce_fg_shared))

colnames(rowData(sce_og_shared))[
  colnames(rowData(sce_og_shared)) == "ID"] <- "ID_og"
colnames(rowData(sce_fg_shared))[
  colnames(rowData(sce_fg_shared)) == "ID"] <- "ID_fg"

colData(sce_og_shared) <- colData(sce_og_shared)[,-20]

# combine
sce_tog <- BiocGenerics::cbind(sce_og_shared, sce_fg_shared)
print(sce_tog)

#-------------------------------------------------------------------------------

# normalize together for visual comparison and marker gene detection
quick_clust <- scran::quickCluster(sce_tog)
sce_tog <- scran::computeSumFactors(sce_tog, cluster = quick_clust)
sce_tog <- scuttle::logNormCounts(sce_tog) 

#-------------------------------------------------------------------------------
# get markergenes = genes that are diffrentially "expressed" (=mapped)

# convert to seurat object
seurat <- Seurat::as.Seurat(sce_tog, counts = "counts", data = "logcounts") 
Seurat::Idents(seurat) <- sce_tog$Genome

# get hvgs for marker genes
gene_var <- scran::modelGeneVar(sce_tog)
hvgs <- scran::getTopHVGs(gene_var, n=nr_hvgs)

# find marker genes for each genome
genomes <- base::unique(sce_tog$Genome)
clustlist <- as.list(genomes)
cluster_markers <- lapply(clustlist, function(x){
  markers <- Seurat::FindMarkers(seurat, 
                                 est.use = "wilcox", 
                                 ident.1 = x, 
                                 features = hvgs, 
                                 logfc.threshold = logFC_sample_dmgs,
                                 min.pct = logFC_sample_dmgs)
  markers <- markers[base::order(base::abs(markers$avg_log2FC),
                                 decreasing=TRUE),]
  markers$which_cluster <- base::rep(x, nrow(markers))
  return(markers)
})

names(cluster_markers) <- genomes
markergenes_og <- cluster_markers[["OneGenome"]]

#-------------------------------------------------------------------------------
base::saveRDS(markergenes_og, snakemake@output[["dmgs"]])

utils::sessionInfo()