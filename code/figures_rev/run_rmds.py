
# output for htmls
# doesn't have to be base
OUTPUT_PATH = "/omics/odcf/analysis/OE0538_projects/DO-0008/data/manuscript1_rev/htmls"

# base path to pass into the .Rmds for saving the pdfs
BASE_PATH = "/omics/odcf/analysis/OE0538_projects/DO-0008/data"

# input and output paths of data and separate plots are separarely defined
# in each .Rmd script for better overview within the script

# all plots geneated using this script are just parts of figures and details
# (but not data) are manually adjusted using Affinity.

#-------------------------------------------------------------------------------

targets = []

targets = targets + [OUTPUT_PATH + "/figure_1.html"]
targets = targets + [OUTPUT_PATH + "/figure_2.html"]
targets = targets + [OUTPUT_PATH + "/figure_3.html"]
targets = targets + [OUTPUT_PATH + "/figure_4.html"]
targets = targets + [OUTPUT_PATH + "/figure_5.html"]

targets = targets + [OUTPUT_PATH + "/figure_s2.html"]
targets = targets + [OUTPUT_PATH + "/figure_s2_alluvial_hsc.html"]
targets = targets + [OUTPUT_PATH + "/figure_s2_alluvial_str.html"]
targets = targets + [OUTPUT_PATH + "/figure_s4.html"]
targets = targets + [OUTPUT_PATH + "/figure_s5.html"]

"""
targets = targets + [OUTPUT_PATH + "/figure_s3.html"]

"""

#-------------------------------------------------------------------------------

localrules: all

rule all: 
  input:
      targets

#-------------------------------------------------------------------------------

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
        mem_mb=40000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_3.html"
    params:
        base_path = BASE_PATH
    threads: 4
    script:
        "figure_3.Rmd"

rule run_fig_4:
    resources:
        mem_mb=80000,
        queues="medium"
    output:
        OUTPUT_PATH + "/figure_4.html"
    params:
        base_path = BASE_PATH
    threads: 4
    conda:
        "../envs/ggpattern.yml"
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
