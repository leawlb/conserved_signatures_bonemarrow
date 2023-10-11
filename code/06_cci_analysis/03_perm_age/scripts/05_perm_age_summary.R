
library(tidyverse)

set.seed(37)

species <- snakemake@params[["species"]]

print(species)

# load res_lists of each species
res_lists <- list()
for(i in 1:length(species)){
  print(species[i])
  print(snakemake@input[["res_list_input_path"]][i])
  if(grepl(species[i], snakemake@input[["res_list_input_path"]][i])){
    res_lists[[i]] <- readRDS(file = snakemake@input[["res_list_input_path"]][i])
    names(res_lists)[i] <- species[i]
  }
}
print(names(res_lists))

large_df_list <- lapply(as.list(species), function(s){
  
  res_list <- res_lists[[s]]
  res_df_list <- list()
  for(i in names(res_list)){
    res_df_list[[i]] <- res_list[[i]]$res_df
    
    res_df_list[[i]]$species <- s
    res_df_list[[i]]$ctp <- i
    res_df_list[[i]]$emitter <- str_split(i, "&")[[1]][1]
    res_df_list[[i]]$receiver <- str_split(i, "&")[[1]][2]
    
    res_df_list[[i]] <- rownames_to_column(res_df_list[[i]], var = "int")
  }
  
  res_df_large <- bind_rows(res_df_list)
  return(res_df_large)
})

summary_df <- bind_rows(large_df_list)

# add more info

# whether an interaction was detected, and is shared or only detected in one age
summary_df$detected <- vector(length = nrow(summary_df))
summary_df$detected[!is.na(summary_df$score_yng) & !is.na(summary_df$score_old)] <- "shared"
summary_df$detected[!is.na(summary_df$score_yng) & is.na(summary_df$score_old)] <- "young only"
summary_df$detected[is.na(summary_df$score_yng) & !is.na(summary_df$score_old)] <- "old only"
summary_df$detected[is.na(summary_df$score_yng) & is.na(summary_df$score_old)] <- "not detected"

print(table(summary_df$detected))

summary_df$detected <- factor(summary_df$detected, 
                             levels = c("not detected", "young only", "shared", "old only"))

# whether a shared interaction delta score is significantly higher or lower than expected from permutation
summary_df$direction <- vector(length = nrow(summary_df))
summary_df$direction[summary_df$p_adj <= 0.05 & summary_df$score_delta > summary_df$perm_scores_delta_median] <- "delta_score significantly higher"
summary_df$direction[summary_df$p_adj <= 0.05 & summary_df$score_delta < summary_df$perm_scores_delta_median] <- "delta_score significantly lower"
summary_df$direction[summary_df$p_adj > 0.05] <- "not significant"
summary_df$direction[is.na(summary_df$p_adj)] <- "pval not estimated"
summary_df$direction[summary_df$detected == "young only"] <- "not shared"
summary_df$direction[summary_df$detected == "old only"] <- "not shared"
summary_df$direction[summary_df$detected == "not detected"] <- "not detected"

print(table(summary_df$direction))

# change with age
summary_df$direction <- factor(summary_df$direction, 
                              levels = c("not detected", "not shared", "delta_score significantly lower", "delta_score significantly higher", "not significant", "pval not estimated"))

summary_df$change_age <- vector(length = nrow(summary_df))
summary_df$change_age[summary_df$direction == "delta_score significantly higher"] <- "increased"
summary_df$change_age[summary_df$direction == "delta_score significantly lower"] <- "decreased"
summary_df$change_age[summary_df$direction == "not significant"] <- "no significant change (detectable)"
summary_df$change_age[summary_df$direction == "pval not estimated"] <- "no significant change (detectable)"
summary_df$change_age[summary_df$detected == "not detected"] <- "not detected"

print(table(summary_df$change_age))

summary_df$change_age <- factor(summary_df$change_age, 
                               levels = c("not detected", "decreased", "increased","no significant change (detectable)"))

# ligands and receptors
summary_df$ligand <- vector(length = nrow(summary_df))
summary_df$receptor <- vector(length = nrow(summary_df))
for(i in 1:nrow(summary_df)){
  summary_df$ligand[i] <- str_split(summary_df$int[i], "_")[[1]][1]
  summary_df$receptor[i] <- str_split(summary_df$int[i], "_")[[1]][2]
}

# remaining factorization
summary_df$species <- factor(summary_df$species, 
                             levels = c("mmus", "mcas", "mspr", "mcar"))

print(table(summary_df$detected, summary_df$shared))

saveRDS(summary_df, snakemake@output[["summary_df"]])