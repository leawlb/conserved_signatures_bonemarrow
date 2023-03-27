
library(DESeq2)
library(tidyverse)
library(ashr)

# in this script, results from DESeq objects are exported for every 
# combination of binary species comparisons
# for differentially expressed genes
#-------------------------------------------------------------------------------

tdsq_list <- readRDS(file = snakemake@input[["tdsq"]])

fc_cutoff <- snakemake@params[["fc_cutoff"]]
padj_cutoff <- snakemake@params[["padj_cutoff"]]
print(fc_cutoff)
print(padj_cutoff)

fraction_curr <- snakemake@wildcards[["fraction"]]
print(fraction_curr)

if(fraction_curr == "hsc"){
  nr_cluster <- 13
}else if(fraction_curr == "str"){
  nr_cluster <- 17
}
print(nr_cluster)
clusters <- paste0("cluster_", c(1:nr_cluster))
print(clusters)

species <- snakemake@params[["species"]]
print(species)

print(snakemake@output[["species_specific_csv"]])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

a <- species[1]
b <- species[2]
c <- species[3]
d <- species[4]

# every combination for easier subsetting later
comb_list<- list(
  c(a, b),
  c(a, c),
  c(a, d),
  c(b, a),
  c(b, c),
  c(b, d),
  c(c, a),
  c(c, b),
  c(c, d),
  c(d, a),
  c(d, b),
  c(d, c)
)

#-------------------------------------------------------------------------------
# RESULTS
# get results for each combination

res_list_list <- lapply(tdsq_list, function(tdsq){
  
  res_list <- lapply(comb_list, function(comb){
    res_lfcs <- lfcShrink(tdsq, contrast=c("condition", comb[[1]], comb[[2]]), 
                          type = "ashr", quiet = TRUE, lfcThreshold = fc_cutoff,
                          alpha = padj_cutoff)
    return(res_lfcs)
  })
  for(i in 1:length(comb_list)){
    names(res_list)[[i]] <- paste0(comb_list[[i]][1], "-", comb_list[[i]][2])
  }
  return(res_list)
})
names(res_list_list) <- c(1:nr_cluster)

#-------------------------------------------------------------------------------
# CATEGORY
# add category information to res dataframes

res_temp1 <- lapply(res_list_list, function(res_list){
  
  res_list_n <- lapply(res_list, function(res){
    
    # remove NAs in log2FC and padj
    res <- res[which(!is.na(res$padj)),]
    res <- res[which(!is.na(res$log2FoldChange)),]
    
    # get category positions
    res$category <- rep("None", n = nrow(res))
    fc_pos <- which(abs(res$log2FoldChange) >= fc_cutoff & res$padj > padj_cutoff)
    pv_pos <- which(abs(res$log2FoldChange) < fc_cutoff & res$padj <= padj_cutoff)
    fp_pos <- which(abs(res$log2FoldChange) >= fc_cutoff & res$padj <= padj_cutoff)
    
    # in rare cases a category doesn't exist
    if(!(length(fc_pos) & length(pv_pos) & length(fp_pos) != 0)){
      print(paste0("fc_pos is ", length(fc_pos)))
      print(paste0("pv_pos is ", length(pv_pos)))
      print(paste0("fp_pos is ", length(fp_pos)))
    }
    
    # add category
    if(length(fc_pos) > 0){res[fc_pos,]$category <- "FC"}
    if(length(pv_pos) > 0){res[pv_pos,]$category <- "pval"}
    if(length(fp_pos) > 0){res[fp_pos,]$category <- "FC + pval"}
    return(res)
  })
  return(res_list_n)
})
names(res_temp1) <- c(1:nr_cluster)

#-------------------------------------------------------------------------------
# REARRANGE 
# rearrange by species for species marker csv export and species-specificity

res_temp2 <- unlist(res_temp1)
names(res_temp2) <- gsub("\\.", "_", names(res_temp2))

# rearrange to species
res_temp3 <- lapply(species, function(s){
  list <- res_temp2[c(grep(s, substr(names(res_temp2), 3, 6)), 
                      grep(s, substr(names(res_temp2), 4, 7)))]
  return(list)
})
names(res_temp3) <- unlist(species)
length(res_temp3[[1]])
print(nr_cluster*3)

# rearrange by cluster within species
res_temp4 <- lapply(res_temp3, function(res_list){
  nr_cluster
  res_list_n <- lapply(c(1:nr_cluster), function(i){
    res <- res_list[grep(paste0("^", i, "[[:punct:]]+"), names(res_list))]
    return(res)
  })
  names(res_list_n) <- c(1:nr_cluster)
  return(res_list_n)
})

#-------------------------------------------------------------------------------
# SPECIFICITY
# add information on species-specificity for each species and cluster

res_temp5 <- lapply(res_temp4, function(res_list_clust){
  
  res_df_clust_n <- lapply(res_list_clust, function(res_list){

    # rearrange dataframe
    for(j in c(1:length(res_list))){
      res_list[[j]] <- as.data.frame(res_list[[j]])
      res_list[[j]]$comparison <- vector(length = nrow(res_list[[j]]))
      res_list[[j]]$comparison <- str_sub(names(res_list)[j], -9,
                                              nchar(names(res_list)[j]))
      res_list[[j]] <- rownames_to_column(res_list[[j]], var = 'gene')
    }
    
    res_df <- rbind(res_list[[1]], res_list[[2]], res_list[[3]])
    
    res_df$specificity <- vector(length = nrow(res_df))
    species_curr <- str_sub(res_df$comparison[1], 1, 4)
    for(i in unique(res_df$gene)){
      if(length(which(res_df$gene == i)) > 3){
        stop("Multiple grepping of genes")
      }
      if(length(which(res_df$category[res_df$gene == i] == "FC + pval")) == 3){
        # if species has significant FC to all other species (positive)
        if(length(which(res_df$log2FoldChange[res_df$gene == i] >= fc_cutoff)) == 3){
          res_df$specificity[res_df$gene == i] <- paste0(species_curr, " only")
          # if species has significant FC to all other species (negative)
        }else if(length(which(res_df$log2FoldChange[res_df$gene == i] < fc_cutoff)) == 3){
          res_df$specificity[res_df$gene == i] <- paste0("all but ", species_curr)
        }
      }
    }    
    
    res_df$species <- rep(species_curr, nrow(res_df))
    
    res_df <- res_df[order(abs(as.numeric(res_df$log2FoldChange)), 
                           decreasing = TRUE),]
    return(res_df)
  })
  
  print(head(res_df_clust_n[[1]]))
  names(res_df_clust_n) <- paste0("cluster_", names(res_df_clust_n))
  return(res_df_clust_n)
})
saveRDS(res_temp5, snakemake@output[["species_res_df"]])

#-------------------------------------------------------------------------------
# make and export list of genes that are specifically (not) expressed in species

lapply(res_temp5, function(res_list_clust){
  
  # get and sort genes
  gene_list <- lapply(res_list_clust, function(cluster_df){
    genes <- sort(unique(cluster_df$gene[cluster_df$specificity != FALSE]))
    print(length(genes))
    return(genes)
  })
  
  nr <- 1000 # number of rows required
  names(gene_list) <- names(res_list_clust)
  build_df <- data.frame(row.names = c(1:nr))
  
  # add genes and required number of empty rows
  for(i in 1:length(gene_list)){
    build_df[,i] <- c(gene_list[[i]], rep("", (nr-length(gene_list[[i]]))))
    colnames(build_df)[i] <- names(gene_list)[i]
  }
  
  species_curr <- res_list_clust[[1]]$species[1]

  # save
  combination <- paste0(fraction_curr, "_", species_curr)
  path_comb <- snakemake@output[["species_specific_csv"]][
    grep(combination, snakemake@output[["species_specific_csv"]])]
  print(path_comb)
  
  write.csv2(build_df, path_comb, row.names = FALSE)
})

#-------------------------------------------------------------------------------
# REARRANGE 
# rearrange back to cluster 

res_temp6 <- list()
for(i in clusters){
  for(j in species){
    res_temp6[[i]][[j]] <- res_temp5[[j]][[i]]
    names(res_temp6)[which(clusters == i)] <- i
    names(res_temp6[[i]])[which(species == j)] <- j 
  }
}
print(names(res_temp6))

#-------------------------------------------------------------------------------
# SUMMARIZE 
# put all entries for one species into one dataframe for all

res_temp7 <- lapply(res_temp6, function(res_list){
  
  res_df <- res_list[[1]]
  for(i in 2:length(res_list)){
    res_df <- rbind(res_df, res_list[[i]])
  }
  print(dim(res_df))
  return(res_df)
})

saveRDS(res_temp7, snakemake@output[["cluster_res_df"]])
