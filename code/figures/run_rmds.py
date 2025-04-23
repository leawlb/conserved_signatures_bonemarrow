
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

targets = targets + [OUTPUT_PATH + "/figures2.html"]


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
        
        
rule run_figs2:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figures2.html"
    threads: 4
    conda:
        "../envs/ggalluvial.yml"
    script:
        "figures1.Rmd"
