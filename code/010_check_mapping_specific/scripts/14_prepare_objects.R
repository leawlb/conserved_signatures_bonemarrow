#-------------------------------------------------------------------------------

library(biomaRt)
library(SingleCellExperiment) 
library(scran) 
library(scater) 
set.seed(37)

#-------------------------------------------------------------------------------

sce_og <- readRDS(file = snakemake@input[["sce_07"]])
sce_fg_path <- snakemake@input[["sce_fg_path"]]

nr_hvgs <- snakemake@params[["nr_hvgs"]]

spec_curr <- snakemake@wildcards[["species"]]
fraction_curr <- snakemake@wildcards[["fraction"]]
individuals <- snakemake@params[["individuals"]]

#-------------------------------------------------------------------------------
# load objects

sce_fg_list <- list()
for(i in individuals){
  if(grepl(spec_curr, i) & grepl(fraction_curr, i)){
    sce_fg_list[[i]] <- readRDS(file = paste0(sce_fg_path, "/", 
                                              spec_curr, "/sce_", i, "-01"))
    # paste individual to colnames for merging later
    colnames(sce_fg_list[[i]]) <- paste0(colnames(sce_fg_list[[i]]), "_", i)
  }
}

#-------------------------------------------------------------------------------
# prepare four genomes objects

# subset barcodes
sce_fg_list_qc <- lapply(sce_fg_list, function(sce_fg){
  sce_fg <- sce_fg[,which(!is.na(match(colnames(sce_fg), colnames(sce_og))))]
  return(sce_fg)
})

# four genomes mspr objects don't have symbols, take from ensembl
if(grepl("mspr", spec_curr)){
  
  ensembl_mspr <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                          host="https://www.ensembl.org",
                          dataset="mspretus_gene_ensembl")
  ensembl_list_mspr=getBM(attributes=c("ensembl_gene_id",
                                       "mmusculus_homolog_ensembl_gene",
                                       "mmusculus_homolog_associated_gene_name"), 
                          mart=ensembl_mspr)
  
  sce_fg_list_qc <- lapply(sce_fg_list_qc, function(sce_fg){
    
    intersect_IDs <- intersect(rowData(sce_fg)$Symbol, 
                               ensembl_list_mspr$ensembl_gene_id)
    
    ensembl_list_mspr <- ensembl_list_mspr[
      ensembl_list_mspr$ensembl_gene_id %in% intersect_IDs,]
    sce_fg <- sce_fg[rowData(sce_fg)$Symbol %in% intersect_IDs,]
    ensembl_list_mspr <- ensembl_list_mspr[
      !duplicated(ensembl_list_mspr$ensembl_gene_id),]
    sce_fg <- sce_fg[!duplicated(rowData(sce_fg)$Symbol),]
    
    rowData(sce_fg)$Symbol[
      match(ensembl_list_mspr$ensembl_gene_id,
            rownames(sce_fg))] <- ensembl_list_mspr$mmusculus_homolog_associated_gene_name
    
    rowData(sce_fg)$Symbol[rowData(sce_fg)$Symbol == ""] <- rowData(sce_fg)$ID[
      rowData(sce_fg)$Symbol == ""]
    
    return(sce_fg)
  })
}

#-------------------------------------------------------------------------------
# prepare both objects for merging
sce_fg_list_shared <- lapply(sce_fg_list_qc, function(sce_fg){
  
  shared_genes <- intersect(rowData(sce_og)$Symbol, rowData(sce_fg)$Symbol)
  length(shared_genes)
  
  sce_fg_shared <- sce_fg[rowData(sce_fg)$Symbol %in% shared_genes,]
  sce_fg_shared <- sce_fg_shared[order(rowData(sce_fg_shared)$Symbol),]
  sce_fg_shared <- sce_fg_shared[!duplicated(rowData(sce_fg_shared)$Symbol),]
  return(sce_fg_shared)
  
})

for(i in 1:length(sce_fg_list_shared)){
  if(i == 1){
    sce_fg_merged <- sce_fg_list_shared[[i]] # set first item of list
  }else{
    sce_fg_merged <- cbind(sce_fg_merged, sce_fg_list_shared[[i]]) 
  }
}


shared_genes <- intersect(rowData(sce_og)$Symbol,
                          rowData(sce_fg_list_shared[[1]])$Symbol)
sce_og_shared <- sce_og[rowData(sce_og)$Symbol %in% shared_genes,]
sce_og_shared <- sce_og_shared[order(rowData(sce_og_shared)$Symbol),]
sce_og_shared <- sce_og_shared[!duplicated(rowData(sce_og_shared)$Symbol),]

sce_fg_merged
sce_og_shared

colnames(sce_og_shared) <- paste0(colnames(sce_og_shared), "_og")
sce_og_shared$Genome <- rep("OneGenome", ncol(sce_og_shared))

colnames(sce_fg_merged) <- paste0(colnames(sce_fg_merged), "_fg")
sce_fg_merged$Genome <- rep("FourGenomes", ncol(sce_fg_merged))

rownames(sce_og_shared) <- rowData(sce_og_shared)$Symbol 
rownames(sce_fg_merged) <- rowData(sce_fg_merged)$Symbol 

colData(sce_og_shared) <- colData(sce_og_shared)[,colnames(colData(sce_og_shared)) %in% colnames(colData(sce_fg_merged))]

colnames(rowData(sce_og_shared))[
  colnames(rowData(sce_og_shared)) == "ID"] <- "ID_og"
colnames(rowData(sce_fg_merged))[
  colnames(rowData(sce_fg_merged)) == "ID"] <- "ID_fg"

assays(sce_og_shared) <- list(counts = counts(sce_og_shared))
assays(sce_fg_merged)$counts[1:10, 1:10]
assays(sce_og_shared)$counts[1:10, 1:10]

reducedDim(sce_og_shared) <- NULL
reducedDim(sce_og_shared) <- NULL

#-------------------------------------------------------------------------------
# merge 
sce_tog <- cbind(sce_og_shared, sce_fg_merged)

# normalize
set.seed(37)

quick_clust <- quickCluster(sce_tog) 
sce_tog <- computeSumFactors(sce_tog, cluster = quick_clust)
sce_tog <- logNormCounts(sce_tog) # log-transformation

sce_tog

saveRDS(sce_tog, file = snakemake@output[["sce_tog"]])
