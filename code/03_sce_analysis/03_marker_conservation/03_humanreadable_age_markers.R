library(tidyverse)


load(
  "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/03_marker_conservation/cons_markers_hsc.RData"
)
markers_cons_hsc <- markers_conservation_hsc
names(markers_cons_hsc)

# load the markers obtained from different ages
load(
  "/omics/odcf/analysis/OE0538_projects/DO-0008/data/scRNAseq/main_analysis/sce_objects/03_sce_analysis/03_marker_conservation/02_age/cons_markers_hsc.RData"
)
markers_conservation_hsc <- markers_conservation_hsc
names(markers_conservation_hsc)


gene_flags <- map_dfr(names(markers_cons_hsc), function(ct) {
  young <- markers_conservation_hsc[[grep(
    paste(ct, "young", sep = "_"),
    names(markers_conservation_hsc),
    value = TRUE
  )]] %>%
    as.data.frame()
  young$gene <- rownames(young)
  young <- young[complete.cases(young), , drop = FALSE]

  old <- markers_conservation_hsc[[grep(
    paste(ct, "old", sep = "_"),
    names(markers_conservation_hsc),
    value = TRUE
  )]] %>%
    as.data.frame()
  old$gene <- rownames(old)
  old <- old[complete.cases(old), , drop = FALSE]

  mark <- markers_cons_hsc[[grep(ct, names(markers_cons_hsc), value = TRUE)]]
  mark <- mark[complete.cases(mark), , drop = FALSE]

  old_not_young_genes <- setdiff(old$gene, young$gene)
  old_not_mark_not_young_genes <- setdiff(
    old$gene,
    union(young$gene, rownames(mark))
  )

  all_genes <- union(old_not_young_genes, old_not_mark_not_young_genes)

  tibble(
    celltype = ct,
    gene = all_genes,
    old_not_young = gene %in% old_not_young_genes,
    old_not_mark_not_young = gene %in% old_not_mark_not_young_genes
  )
})

View(gene_flags)





