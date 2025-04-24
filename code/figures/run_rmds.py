
# output for htmls
# doesn't have to be base
OUTPUT_PATH = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/manuscript1/htmls"

# input and output paths of data and separate plots are separarely defined
# in each .Rmd script for better overview within the script

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_PATH + "/figure1.html"]
targets = targets + [OUTPUT_PATH + "/figure2.html"]
targets = targets + [OUTPUT_PATH + "/figure3.html"]
targets = targets + [OUTPUT_PATH + "/figure4.html"]
targets = targets + [OUTPUT_PATH + "/figure5.html"]
targets = targets + [OUTPUT_PATH + "/figure6.html"]

targets = targets + [OUTPUT_PATH + "/figure_s2.html"]
targets = targets + [OUTPUT_PATH + "/figure_s3.html"]
targets = targets + [OUTPUT_PATH + "/figure_s4.html"]
#targets = targets + [OUTPUT_PATH + "/figure_s3_silhouette.html"]
#targets = targets + [OUTPUT_PATH + "/figure_s3_genes_v_pseudotime.html"]


#-------------------------------------------------------------------------------

localrules: all  

rule all: 
  input:
      targets

rule run_fig1:
    resources:
        mem_mb=80000,
        queues = "medium-debian"
    output:
        OUTPUT_PATH + "/figure1.html",
    threads: 4
    script:
        "figure1.Rmd"
        
rule run_fig2:
    resources:
        mem_mb=80000,
        queues = "medium-debian"
    output:
        OUTPUT_PATH + "/figure2.html",
    threads: 4
    conda:
        "../envs/ggpattern.yml"
    script:
        "figure2.Rmd"
  
rule run_fig3:
    resources:
        mem_mb=80000,
        queues = "medium-debian"
    output:
        OUTPUT_PATH + "/figure3.html",
    threads: 4
    conda:
        "../envs/ggpattern.yml"
    script:
        "figure3.Rmd"

rule run_fig4:
    resources:
        mem_mb=80000,
        queues = "medium-debian"
    output:
        OUTPUT_PATH + "/figure4.html",
    threads: 4
    script:
        "figure4.Rmd"

    
rule run_fig5:
    resources:
        mem_mb=80000,
        queues = "medium-debian"
    output:
        OUTPUT_PATH + "/figure5.html",
    threads: 4
    script:
        "figure5.Rmd"

rule run_fig6:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure6.html"
    threads: 4
    conda:
        "../envs/upsetplot.yml"
    script:
        "figure6.Rmd"
        
        
rule run_fig_s2:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure_s2.html"
    threads: 4
    conda:
        "../envs/ggalluvial.yml"
    script:
        "figure_s2.Rmd"
        
rule run_fig_s3:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure_s3.html"
    threads: 4
    script:
        "figure_s3.Rmd"
        
rule run_fig_s3_silhouette:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure_s3_silhouette.html"
    threads: 4
    script:
        "figure_s3_silhouette.Rmd"
        
rule run_fig_s3_genes_v_pseudotime:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure_s3_genes_v_pseudotime.html"
    threads: 4
    script:
        "figure_s3_genes_v_pseudotime.Rmd"

rule run_fig_s4:
    resources:
        mem_mb=200000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure_s4.html"
    threads: 4
    script:
        "figure_s4.Rmd"
