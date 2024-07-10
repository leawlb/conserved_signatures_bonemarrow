#-------------------------------------------------------------------------------
# unify list of all dmgs to exclude from SCE objects

dmg_paths <- snakemake@input[["dmgs"]]

dmg_list <- list()
for(i in 1:length(dmg_paths)){
  dmg_list[[i]] <- base::readRDS(file = dmg_paths[[i]])
}
print(dmg_list)

dmgs_exclude <- rownames(dmg_list[[1]])

for(i in 2:length(dmg_paths)){
  dmgs_exclude <- c(dmgs_exclude, rownames(dmg_list[[i]]))
}

dmgs_exclude <- base::unique(dmgs_exclude)
dmgs_exclude <- base::sort(dmgs_exclude)
print(dmgs_exclude)

base::saveRDS(dmgs_exclude, file = snakemake@output[["dmg_list"]])

utils::sessionInfo()