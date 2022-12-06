
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mcar', '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mcas', '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mmus', '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mspr', "sce_04" = c('/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mcar', '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mcas', '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mmus', '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/04_norm/mspr')),
    output = list('/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/06_mrge/sce_all-06', "sce_06" = '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects/06_mrge/sce_all-06'),
    params = list(c('mcar_yng_hsc_1_0', 'mcar_old_str_1_0', 'mcar_old_str_2_0', 'mcar_old_hsc_3_0', 'mcar_old_hsc_1_0', 'mcar_old_hsc_2_0', 'mcar_old_str_3_1', 'mcar_old_str_3_2', 'mcar_old_str_5_1', 'mcar_yng_str_1_0', 'mcar_old_str_5_2', 'mcar_old_hsc_4_0', 'mcar_old_hsc_5_0', 'mcar_yng_hsc_3_0', 'mcar_yng_str_2_2', 'mcar_yng_str_2_1', 'mcar_yng_hsc_2_0', 'mcar_yng_str_2_3', 'mcas_old_hsc_1_0', 'mcas_old_hsc_3_0', 'mcas_old_str_1_0', 'mcas_old_hsc_2_0', 'mcas_old_str_2_0', 'mcas_old_str_3_0', 'mcas_yng_hsc_1_0', 'mcas_yng_str_1_0', 'mcas_yng_str_2_0', 'mcas_yng_hsc_2_0', 'mcas_yng_hsc_3_0', 'mcas_yng_str_3_0', 'mcas_yng_hsc_4_0', 'mcas_old_hsc_5_0', 'mcas_old_hsc_4_0', 'mcas_old_str_4_0', 'mcas_old_str_5_0', 'mcas_yng_str_6_0', 'mcas_yng_hsc_5_0', 'mcas_yng_hsc_6_0', 'mmus_old_hsc_3_0', 'mmus_old_hsc_1_0', 'mmus_old_str_2_0', 'mmus_old_str_3_0', 'mmus_old_str_1_0', 'mmus_old_hsc_2_0', 'mmus_yng_hsc_1_0', 'mmus_yng_str_1_0', 'mmus_yng_hsc_2_0', 'mmus_yng_str_2_0', 'mmus_yng_str_3_0', 'mmus_yng_hsc_3_0', 'mmus_old_hsc_4_0', 'mspr_old_str_2_0', 'mspr_old_hsc_1_0', 'mspr_old_str_1_0', 'mspr_old_hsc_2_0', 'mspr_old_hsc_3_0', 'mspr_old_str_3_0', 'mspr_old_hsc_4_0', 'mspr_old_str_5_0', 'mspr_old_str_4_0', 'mspr_old_hsc_5_0', 'mspr_yng_str_1_0', 'mspr_yng_hsc_2_0', 'mspr_yng_str_3_2', 'mspr_yng_hsc_3_0', 'mspr_yng_str_3_1'), c('mcar', 'mcas', 'mmus', 'mspr'), "individuals" = c('mcar_yng_hsc_1_0', 'mcar_old_str_1_0', 'mcar_old_str_2_0', 'mcar_old_hsc_3_0', 'mcar_old_hsc_1_0', 'mcar_old_hsc_2_0', 'mcar_old_str_3_1', 'mcar_old_str_3_2', 'mcar_old_str_5_1', 'mcar_yng_str_1_0', 'mcar_old_str_5_2', 'mcar_old_hsc_4_0', 'mcar_old_hsc_5_0', 'mcar_yng_hsc_3_0', 'mcar_yng_str_2_2', 'mcar_yng_str_2_1', 'mcar_yng_hsc_2_0', 'mcar_yng_str_2_3', 'mcas_old_hsc_1_0', 'mcas_old_hsc_3_0', 'mcas_old_str_1_0', 'mcas_old_hsc_2_0', 'mcas_old_str_2_0', 'mcas_old_str_3_0', 'mcas_yng_hsc_1_0', 'mcas_yng_str_1_0', 'mcas_yng_str_2_0', 'mcas_yng_hsc_2_0', 'mcas_yng_hsc_3_0', 'mcas_yng_str_3_0', 'mcas_yng_hsc_4_0', 'mcas_old_hsc_5_0', 'mcas_old_hsc_4_0', 'mcas_old_str_4_0', 'mcas_old_str_5_0', 'mcas_yng_str_6_0', 'mcas_yng_hsc_5_0', 'mcas_yng_hsc_6_0', 'mmus_old_hsc_3_0', 'mmus_old_hsc_1_0', 'mmus_old_str_2_0', 'mmus_old_str_3_0', 'mmus_old_str_1_0', 'mmus_old_hsc_2_0', 'mmus_yng_hsc_1_0', 'mmus_yng_str_1_0', 'mmus_yng_hsc_2_0', 'mmus_yng_str_2_0', 'mmus_yng_str_3_0', 'mmus_yng_hsc_3_0', 'mmus_old_hsc_4_0', 'mspr_old_str_2_0', 'mspr_old_hsc_1_0', 'mspr_old_str_1_0', 'mspr_old_hsc_2_0', 'mspr_old_hsc_3_0', 'mspr_old_str_3_0', 'mspr_old_hsc_4_0', 'mspr_old_str_5_0', 'mspr_old_str_4_0', 'mspr_old_hsc_5_0', 'mspr_yng_str_1_0', 'mspr_yng_hsc_2_0', 'mspr_yng_str_3_2', 'mspr_yng_hsc_3_0', 'mspr_yng_str_3_1'), "species" = c('mcar', 'mcas', 'mmus', 'mspr')),
    wildcards = list(),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("paths" = list("output_dir" = '/omics/odcf/analysis/OE0538_projects/DO-0008/data/sce_objects', "cellranger_output" = '/01_cellranger_output', "preprocessing" = '/02_preprocessing', "integration" = '/03_integration'), "metadata" = list("raw" = '/omics/odcf/analysis/OE0538_projects/DO-0008/metadata/OE0538_DO-0008_metadata_combined.csv', "identifiers" = c('Species_ID', 'Age_ID', 'Fraction_ID', 'Sample_NR'), "values" = list("cutoff_umis" = 100, "cutoff_doublets" = 7.5, "cutoff_mitos" = 5, "cutoff_sum" = 400, "cutoff_detected" = 300, "nr_hvgs" = 2000, "nr_hvgs_batch_correction" = 5000), "batches" = c('individual', 'Batch_exp_day', 'Batch_sequencing'), "batch_use" = c('individual')), "run_preprocessing_summary" = FALSE, "run_mnncorrect" = FALSE, "run_seurat3" = FALSE, "run_scmerge" = FALSE, "hvgs_for_batch_correction" = FALSE, "rescale_for_batch_correction" = TRUE),
    rule = 'merge_datasets_all',
    bench_iteration = as.numeric(NA),
    scriptdir = '/home/l012t/repositories/Interspecies_BM_phd/code/03_integration/scripts',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
#-------------------------------------------------------------------------------

library(stringr)

sce_04_path <- snakemake@input[["sce_04"]] # correct preprocessing output
sce_06_path <- snakemake@output[["sce_06"]] # correct integration output
individuals <- snakemake@params[["individuals"]] # all objects to merge

species_curr <- snakemake@wildcards[["species"]]

sce_list_04 <- list()

if(!is.null(species_curr)){ # if species wildcard was used (merge_species)
  for(i in individuals){
    if(grepl(species_curr, i) == TRUE){
      sce_list_04[[i]] <- readRDS(file = paste0(sce_04_path, "/sce_", i, "-04"))
    }
  }
}

if(is.null(species_curr)){ # if no wildcard was used (merge all)
  for(j in 1:length(sce_04_path)){
    # the "current" species i is the last part of the path
    species_curr <- str_sub(sce_04_path[j], start= -4) 
    for(i in individuals){
      if(grepl(species_curr, i) == TRUE){
        sce_list_04[[i]] <- readRDS(file = paste0(sce_04_path[j], 
                                                  "/sce_", i, "-04"))
      }
    }
  }
}

print(names(sce_list_04))

#-------------------------------------------------------------------------------
## prepare merge

# get genes shared by all objects
rownames_list <- lapply(sce_list_04, function(sce){
  
  print(rownames(sce)[1]) # this only works with this line of code here? why? DO NOT DELETE 
  rn <- rownames(sce)
  return(rn)
})
shared_genes <- Reduce(intersect, rownames_list)

# subset all objects by shared_genes
sce_list_04 <- lapply(sce_list_04, function(sce){
  sce <- sce[rownames(sce) %in% shared_genes,]
  return(sce)
})

# add individual sample name to barcodes to avoid random doublets
sce_list_04 <- lapply(sce_list_04, function(sce){
  colnames(sce) <- paste0(colnames(sce), "_", sce$Object_ID)
  return(sce)
})

#-------------------------------------------------------------------------------
## all merge  

# discard designated objects as specified in metadata
# only merge all objects if specified in config
for(i in 1:length(sce_list_04)){
  if(i == 1){
    if(sce_list_04[[i]]$Keep_sample[1] == TRUE){ 
        sce_merged <- sce_list_04[[i]] # set first item of list
      }else{ # only if not to be discarded
        stop("first item in list is to be discarded")
      }
    }else{
      if(sce_list_04[[i]]$Keep_sample[1] == TRUE){
        sce_merged <- cbind(sce_merged, sce_list_04[[i]]) 
    }
  }
}
print(nrow(sce_merged))

saveRDS(sce_merged, file = sce_06_path)