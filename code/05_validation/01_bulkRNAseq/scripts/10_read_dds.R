
library(tidyverse)
library(DESeq2)

# generate a DESeq2 object

#------------------------------------------------------------------------------
# genome annotation
# ensembl for BL6 has been downloaded for scRNaseq data before, see metadata

ensembl_path <- snakemake@input[["ensembl_mus_path"]]
ensembl_mus <- base::readRDS(ensembl_path)

#------------------------------------------------------------------------------
# counts (from FeatureCounts)

counts_path <- snakemake@input[["counts_path"]]

counts <- read_delim(
  counts_path, 
  delim = "\t",
  skip = 1)

colnames(counts)

# change the names of the entire path to the actual sample name
counts <- counts %>%
  as.data.frame() %>%
  dplyr::rename_with(~ str_replace(.x, "^.*/([^.]*)\\..*$", "\\1")) %>% 
  dplyr::rename_with(~ sub("_[^_]*_[^_]*$", "", .x)) %>%
  dplyr::rename_with(~ sub("^[^_]*_", "", .x)) 

# reorder as counts table
counts <- counts[,c(1, 7:ncol(counts))]

counts <- column_to_rownames(counts, var = "Geneid")

#------------------------------------------------------------------------------
# metadata

# this is the DKFZ-style metadata table from sequencing submission
metadata_path <- snakemake@input[["metadata_path"]]

metadata <- read_delim(
  metadata_path, 
  col_types = "c",
  delim = ",")

metadata <- metadata[,
  c(
    "SAMPLE_NAME_GPCF",
    "PID",
    "INDIVIDUAL",
    "CAGE_ID",
    "GENDER",
    "Species Scientific Name",
    "GENOTYPE",
    "DATE_OF_BIRTH",
    "DATE_OF_DEATH",
    "Number of Reads",
    "CELL_INPUT[TOTAL_ALIVE CELLS]"
    )] %>%
  dplyr::distinct()

# rename for consistency 
colnames(metadata) <- c(
    "Object_ID",
    "Patient_ID",
    "Individual",
    "Cage_ID",
    "Sex",
    "Species",
    "Genotype",
    "Date_of_birth",
    "Date_of_death",
    "Read_number",
    "Cell_number")

# order by sample name in counts object
metadata <- metadata[match(
  metadata$Object_ID, 
  colnames(counts)),]

# add additional metadata encoded in the sample name
metadata <- metadata %>% 
  separate(
    Object_ID,
    into = c("Species_ID", "Fraction_ID", "Target_ID", "Condition_ID", "Mouse_ID"),
    sep = "_",
    remove = FALSE
  ) 

#------------------------------------------------------------------------------
# put them all together

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Mouse_ID
  # the design can be changed later
)

#------------------------------------------------------------------------------
# add anno info to rowData
rowData(dds)$gene_id <- ensembl_mus$ensembl_gene_id[
  match(rownames(dds), ensembl_mus$ensembl_gene_id)
]

rowData(dds)$gene_name <- ensembl_mus$external_gene_name[
  match(rownames(dds), ensembl_mus$ensembl_gene_id)
]

rownames(dds) <- rowData(dds)$gene_name

# remove duplicates or empty rows.
dds <- dds[!is.na(rownames(dds)),]

table(duplicated(rownames(dds)))
dds <- dds[!duplicated(rownames(dds)),]

#------------------------------------------------------------------------------
# factoring is required for DESeq2 and also nicer

dds$Object_ID <- factor(
  dds$Object_ID,
  levels = c(
    "mmus_msc_neo1_n_05",
    "mmus_msc_neo1_n_06",
    "mmus_msc_neo1_n_07",
    "mmus_msc_neo1_n_08",
    "mmus_msc_neo1_p_06",
    "mmus_msc_neo1_p_07",
    "mmus_msc_neo1_p_08",
    "mmus_msc_trkb_n_05",
    "mmus_msc_trkb_n_06",
    "mmus_msc_trkb_n_07",
    "mmus_msc_trkb_n_08",
    "mmus_msc_trkb_p_05",
    "mmus_msc_trkb_p_06",
    "mmus_msc_trkb_p_07",
    "mmus_msc_trkb_p_08"
  ))

dds$Target_ID <- factor(
  dds$Target_ID,
  levels = c(
    "neo1",
    "trkb"
  ))

dds$Condition_ID <- factor(
  dds$Condition_ID,
  levels = c(
    "n",
    "p"
  ))

dds$Mouse_ID <- factor(
  dds$Mouse_ID,
  levels = c(
    "05",
    "06",
    "07",
    "08"
  ))

output_path <- snakemake@output[["output_path_dds"]]

saveRDS(dds, output_path)

# check one important gene as example
"Cxcl12" %in% rowData(dds)$gene_name

sessionInfo()