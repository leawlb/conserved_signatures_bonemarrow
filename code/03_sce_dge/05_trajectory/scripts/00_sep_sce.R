
library(scater)

set.seed(37)

sce <- readRDS(snakemake@input[["sce_input"]])
species <- snakemake@params[["species"]]

print(species)
print(snakemake@output[["sce_output_path"]])

for(s in species){
  sce_spc <- sce[,sce$Species_ID == s]
  
  output_path <- snakemake@output[["sce_output_path"]][[grep(s, snakemake@output[["sce_output_path"]])]]
  saveRDS(sce_spc, output_path)
}
