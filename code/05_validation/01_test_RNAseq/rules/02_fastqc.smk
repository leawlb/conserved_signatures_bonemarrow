
# FASTQC for generating a QC html for each fastq file

# intended as an overview of the quality BEFORE any pre-processing
# using this snakemake wrapper:
# https://github.com/s-andrews/FastQC

rule fastqc_raw:
    input:
        # for actual use
        #rules.link_files.output
        # for TESTING!! the file should be "renamed"
        OUTPUT_BASE + "/raw_msc_test/{sample_folder_name}/{sample}_R{read}.fastq.gz"
    output:
        html = OUTPUT_BASE + "/reads/{sample_folder_name}/02_prepro_qc/02_{sample}_R{read}_raw.html",
        zip = OUTPUT_BASE + "/reads/{sample_folder_name}/02_prepro_qc/02_{sample}_R{read}_raw_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet",
        queue = "short" # for DKFZ LSF cluster queue 
    log:
        "logs/{sample_folder_name}/02_fastqc_raw/02_{sample}_R{read}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v5.7.0/bio/fastqc"