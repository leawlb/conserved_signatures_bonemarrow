
# snakemake wrapper 
# FASTQC for generating a QC html for each fastq file
# https://github.com/s-andrews/FastQC

rule fastqc:
    input:
        rules.link_files.output
    output:
        html = OUTPUT_BASE + "/reads/{sample_folder_name}/02_prepro_qc/{sample_r}.html",
        zip = OUTPUT_BASE + "/reads/{sample_folder_name}/02_prepro_qc/{sample_r}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet",
        queue = "medium" # for DKFZ LSF cluster queue 
    log:
        "logs/02_prepro_qc/{sample_folder_name}/{sample_r}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v5.7.0/bio/fastqc"