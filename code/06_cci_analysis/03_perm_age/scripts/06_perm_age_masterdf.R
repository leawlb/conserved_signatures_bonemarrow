
# load perm_age_summary_df
summary_df <- readRDS(file = snakemake@input[["summary_df"]])

# load lots of different data
ident_pair_info_p <- snakemake@input[["ident_pair_info"]]
ident_pair_info_p <- list("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/10_ipis/ipi_mcar_old",
                          "/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/10_ipis/ipi_mmus_yng")
ident_info_p <- snakemake@input[["ident_info"]]
ident_info_p <- list("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/11_idis/idi_mcar_old",
                          "/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/11_idis/idi_mmus_yng")
ident_nrlrs_info_p <- snakemake@input[["ident_nrlrs_info"]]
ident_nrlrs_info_p <- list("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/11_idis/nrlrs_mcar_old",
                          "/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/11_idis/nrlrs_mmus_yng")
sce_p <- snakemake@input[["sce_input"]]
sce_p <- list("/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/01_sprp/sce_mcar_old-01",
                          "/omics/odcf/analysis/OE0538_projects/DO-0008/data/main_analysis/cci_objects/01_cci_preparation/01_sprp/sce_mmus_yng-01")

ident_pair_info_list <- list()
ident_info_list <- list()
ident_nrlrs_info_list <- list()
sce_list <- list()

for(i in 1:length(ident_pair_info_p)){
  ident_pair_info_list[[i]] <- readRDS(file = ident_pair_info_p[[i]])
  ident_info_list[[i]] <- readRDS(file = ident_info_p[[i]])
  ident_nrlrs_info_list[[i]] <- readRDS(file = ident_nrlrs_info_p[[i]])
  sce_list[[i]] <- readRDS(file = sce_p[[i]])
}

ipio_df <- data.frame()
idio_df <- data.frame()
nrlrso_df <- data.frame()
sce <- sce_list[[i]][,0]

# put all into 1 DF
for(i in 1:length(ident_pair_info_list)){
  ipio_df <- rbind(ipio_df, ident_pair_info_list[[i]]$overview)
  idio_df <- rbind(idio_df, ident_info_list[[i]]$overview)
  nrlrso_df <- rbind(nrlrso_df, ident_nrlrs_info_list[[i]])
  sce <- cbind(sce, sce_list[[i]])
}

colnames(colData(sce))

species <- unfactor(unique(summary_df$species))
ctps <- unique(summary_df$ctp)
emitters <- unique(summary_df$emitter)
receivers <- unique(summary_df$receiver)

summary_df$nr_ints_ctp_yng <- vector(length = nrow(summary_df))
summary_df$nr_ints_ctp_old <- vector(length = nrow(summary_df))
for(s in species){
  for(c in ctps){
    summary_df$nr_ints_ctp_yng[summary_df$species == s & summary_df$ctp == c] <- ipio_df$nr_ints[ipio_df$species == s & ipio_df$ident_pairs == c & ipio_df$age == "yng"]
    summary_df$nr_ints_ctp_old[summary_df$species == s & summary_df$ctp == c] <- ipio_df$nr_ints[ipio_df$species == s & ipio_df$ident_pairs == c & ipio_df$age == "old"]
  }
}

summary_df$nr_ints_emi_yng <- vector(length = nrow(summary_df))
summary_df$nr_ints_emi_old <- vector(length = nrow(summary_df))
summary_df$nr_nlrs_emi_yng <- vector(length = nrow(summary_df))
summary_df$nr_nlrs_emi_old <- vector(length = nrow(summary_df))
summary_df$nr_cells_emi_yng <- vector(length = nrow(summary_df))
summary_df$nr_cells_emi_old <- vector(length = nrow(summary_df))
for(s in species){
  for(e in emitters){
    summary_df$nr_ints_emi_yng[summary_df$species == s & summary_df$emitter == e] <- idio_df$nr_ints[idio_df$species == s & idio_df$identity == e & idio_df$age == "yng"]
    summary_df$nr_ints_emi_old[summary_df$species == s & summary_df$emitter == e] <- idio_df$nr_ints[idio_df$species == s & idio_df$identity == e & idio_df$age == "old"]
    
    summary_df$nr_nlrs_emi_yng[summary_df$species == s & summary_df$emitter == e] <- nrlrso_df$nr_lrs[nrlrso_df$species == s & nrlrso_df$identity == e & nrlrso_df$age == "yng"]
    summary_df$nr_nlrs_emi_old[summary_df$species == s & summary_df$emitter == e] <- nrlrso_df$nr_lrs[nrlrso_df$species == s & nrlrso_df$identity == e & nrlrso_df$age == "old"]
    
    summary_df$nr_cells_emi_yng[summary_df$species == s & summary_df$emitter == e] <- length(grep(e, sce$Identity[sce$Species_ID == s & sce$Age_ID == "yng"]))
    summary_df$nr_cells_emi_old[summary_df$species == s & summary_df$emitter == e] <- length(grep(e, sce$Identity[sce$Species_ID == s & sce$Age_ID == "old"]))
  }
}

summary_df$nr_ints_rec_yng <- vector(length = nrow(summary_df))
summary_df$nr_ints_rec_old <- vector(length = nrow(summary_df))
summary_df$nr_nlrs_rec_yng <- vector(length = nrow(summary_df))
summary_df$nr_nlrs_rec_old <- vector(length = nrow(summary_df))
summary_df$nr_cells_rec_yng <- vector(length = nrow(summary_df))
summary_df$nr_cells_rec_old <- vector(length = nrow(summary_df))
for(s in species){
  for(r in receivers){
    summary_df$nr_ints_rec_yng[summary_df$species == s & summary_df$receiver == r] <- idio_df$nr_ints[idio_df$species == s & idio_df$identity == r & idio_df$age == "yng"]
    summary_df$nr_ints_rec_old[summary_df$species == s & summary_df$receiver == r] <- idio_df$nr_ints[idio_df$species == s & idio_df$identity == r & idio_df$age == "old"]
    
    summary_df$nr_nlrs_rec_yng[summary_df$species == s & summary_df$receiver == r] <- nrlrso_df$nr_lrs[nrlrso_df$species == s & nrlrso_df$identity == r & nrlrso_df$age == "yng"]
    summary_df$nr_nlrs_rec_old[summary_df$species == s & summary_df$receiver == r] <- nrlrso_df$nr_lrs[nrlrso_df$species == s & nrlrso_df$identity == r & nrlrso_df$age == "old"]
    
    summary_df$nr_cells_rec_yng[summary_df$species == s & summary_df$receiver == r] <- length(grep(r, sce$Identity[sce$Species_ID == s & sce$Age_ID == "yng"]))
    summary_df$nr_cells_rec_old[summary_df$species == s & summary_df$receiver == r] <- length(grep(r, sce$Identity[sce$Species_ID == s & sce$Age_ID == "old"]))
  }
}

saveRDS(summary_df, file = snakemake@input[["master_df"])
