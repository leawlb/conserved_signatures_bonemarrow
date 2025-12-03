
# output for htmls
# doesn't have to be base
OUTPUT_PATH = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/manuscript1_rev/htmls"

# input and output paths of data and separate plots are separarely defined
# in each .Rmd script for better overview within the script

# all plots geneated using this script are just parts of figures and details
# (but not data) are manually adjusted using Affinity.

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_PATH + "/figure_s3.html"]
targets = targets + [OUTPUT_PATH + "/figure_s4.html"]

#-------------------------------------------------------------------------------

localrules: all

rule all: 
  input:
      targets

rule run_fig_s3:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure_s3.html"
    conda:
        "../envs/ggalluvial.yml"
    threads: 4
    script:
        "figure_s3.Rmd"

rule run_fig_s4:
    resources:
        mem_mb=80000,
        queues="medium-debian"
    output:
        OUTPUT_PATH + "/figure_s4.html"
    conda:
        "../envs/ggalluvial.yml"
    threads: 4
    script:
        "figure_s4.Rmd"