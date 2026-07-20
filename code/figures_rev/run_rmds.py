

# base path to data repository
BASE_PATH = config["base"] + "/data"

# output directory for htmls
OUTPUT_PATH = BASE_PATH + "/manuscript1_rev/htmls"

# input paths for data and 
# output paths of figure pdfs 
# are separarely define in each .Rmd script for better overview of
# where the data came from, and where it goes.

# all plots generated using this script are parts of figures.
# Figures are manually assembled using Affinity (placement, orientation, text).
# Of course, actual data is not changed in Affinity.

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_PATH + "/figure_1.html"]
targets = targets + [OUTPUT_PATH + "/figure_2.html"]
targets = targets + [OUTPUT_PATH + "/figure_3.html"]
targets = targets + [OUTPUT_PATH + "/figure_4.html"]
targets = targets + [OUTPUT_PATH + "/figure_5.html"]
targets = targets + [OUTPUT_PATH + "/figure_6.html"]
# fig 7 is assembled from different other Rmds, see below

# fig s1 is exclusively made from FlowJo and Affinity
targets = targets + [OUTPUT_PATH + "/figure_s2.html"]
targets = targets + [OUTPUT_PATH + "/figure_s2_alluvial_hsc.html"]
targets = targets + [OUTPUT_PATH + "/figure_s2_alluvial_str.html"]
targets = targets + [OUTPUT_PATH + "/figure_s3.html"]
targets = targets + [OUTPUT_PATH + "/figure_s3_genes.html"]
targets = targets + [OUTPUT_PATH + "/figure_s3_silhouette.html"]
targets = targets + [OUTPUT_PATH + "/figure_s3_genes_v_pseudotime.html"]
targets = targets + [OUTPUT_PATH + "/figure_s4.html"]
targets = targets + [OUTPUT_PATH + "/figure_s5.html"]
targets = targets + [OUTPUT_PATH + "/figure_s6.html"]
targets = targets + [OUTPUT_PATH + "/figure_s7_facs.html"]

# other figure parts (thematic)
targets = targets + [OUTPUT_PATH + "/figure_thresholds.html"]
# the rest is run manually for now

#-------------------------------------------------------------------------------

localrules: all

rule all: 
  input:
      targets

#-------------------------------------------------------------------------------
# main figures

rule run_fig_1:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_1.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_1.Rmd"

rule run_fig_2:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_2.html"
    params:
        base_path = BASE_PATH
    threads: 4
    conda:
        "../envs/ggpattern.yml"
    script:
        "figure_2.Rmd"

rule run_fig_3:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_3.html"
    params:
        base_path = BASE_PATH
    threads: 4
    conda:
        "../envs/ggpattern.yml"
    script:
        "figure_3.Rmd"

rule run_fig_4:
    resources:
        mem_mb=80000,
        queues="short"
    output:
        OUTPUT_PATH + "/figure_4.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_4.Rmd"

rule run_fig_5:
    resources:
        mem_mb=80000,
        queues="short"
    output:
        OUTPUT_PATH + "/figure_5.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_5.Rmd"

rule run_fig_6:
    resources:
        mem_mb=80000,
        queues="short"
    output:
        OUTPUT_PATH + "/figure_6.html"
    params:
        base_path = BASE_PATH
    threads: 4
    conda:
        "../envs/bulk_figures5.yaml"
    script:
        "figure_6.Rmd"

#-------------------------------------------------------------------------------
# supplementary figures

rule run_fig_s2:
    resources:
        mem_mb=120000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s2.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_s2.Rmd"

rule run_fig_s2_alluvial_hsc:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s2_alluvial_hsc.html"
    params:
        base_path = BASE_PATH
    conda:
        "../envs/ggalluvial.yml"
    threads: 4
    script:
        "figure_s2_alluvial_hsc.Rmd"

rule run_fig_s2_alluvial_str:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s2_alluvial_str.html"
    params:
        base_path = BASE_PATH
    conda:
        "../envs/ggalluvial.yml"
    threads: 4
    script:
        "figure_s2_alluvial_str.Rmd"

rule run_fig_s3:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s3.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_s3.Rmd"

rule run_fig_s3_genes:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s3_genes.html"
    params:
        base_path = BASE_PATH
    threads: 4
    conda:
        "../envs/ggalluvial.yml" 
        # ggalluvial not required but other tidyverse versions threw errors
    script:
        "figure_s3_genes.Rmd"

rule run_fig_s3_silhouette:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s3_silhouette.html"
    params:
        base_path = BASE_PATH
    threads: 4
    conda:
        "../envs/ggpattern.yml" 
        # ggpatterns not required with these versions seurat conversion works
    script:
        "figure_s3_silhouette.Rmd"

rule run_fig_s3_genes_v_pseudotime:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s3_genes_v_pseudotime.html"
    params:
        base_path = BASE_PATH
    threads: 4
    conda:
        "../envs/ggpattern.yml" 
        # ggpatterns not required with these versions seurat conversion works
    script:
        "figure_s3_genes_v_pseudotime.Rmd"

rule run_fig_s4:
    resources:
        mem_mb=120000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s4.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_s4.Rmd"

rule run_fig_s5:
    resources:
        mem_mb=50000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s5.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_s5.Rmd"

rule run_fig_s6:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_s6.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_s6.Rmd"

rule run_fig_s7:
    resources:
        mem_mb=40000,
        queues="short"
    output:
        OUTPUT_PATH + "/figure_s7_facs.html"
    params:
        base_path = BASE_PATH
    threads: 1
    conda:
        "../envs/bulk_figures5.yaml"
    script:
        "figure_s7_facs.Rmd"

#-------------------------------------------------------------------------------
# there are additional files not sorted by figure, but by topic:

# figure_thresholds
# it was not immediately clear where it would go, but is part of Figure S3 now.
rule run_thresholds:
    resources:
        mem_mb=20000,
        queues="short"
    output:
        OUTPUT_PATH + "/figure_thresholds.html"
    params:
        base_path = BASE_PATH
    threads: 1
    script:
        "figure_thresholds.Rmd"

# motorcortex_mouse_AUCell.R
# .R script, not included in this snakemake pipeline; run manually 
# AND INSERT CORRECT BASE PATH
# data to visualise has been generated in directory 08_sce_brain

# mouse_lemur_AUCell.Rmd
# it was not immediately clear where it would go, but is part of Figure 7 and S8 now.
# run manually 
# AND INSERT CORRECT BASE PATH

# reclaiming_sigs.Rmd
# these plots are part of several figures (figure S4, S5 and S6)
# not all plots were ultimately used
# run manually
# AND INSERT CORRECT BASE PATH
# SOME DATA IS GENERATED AND SAVED HERE: 04_signatures/03_custom_clustering 

# tabula_sapiens_AUCell.Rmd
# this is not currently part of the manuscript
# run manually 
# AND INSERT CORRECT BASE PATH

# zebrafish_AUCell
# it was not immediately clear where it would go, but is part of Figure S8 now.
# run manually 
# AND INSERT CORRECT BASE PATH