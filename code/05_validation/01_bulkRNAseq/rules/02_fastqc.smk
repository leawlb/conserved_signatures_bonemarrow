
# FASTQC for generating a QC html for each fastq file

# intended as an overview of the quality BEFORE any pre-processing
# using this snakemake wrapper:
# https://snakemake-wrappers.readthedocs.io/en/v5.7.0/wrappers/fastqc.html

rule fastqc_raw:
    input:
        # for actual use
        rules.link_files.output
    output:
        html = OUTPUT_BASE + "/reads/{sample_folder_name}/02_fastqc_raw/02_{sample_r}_raw.html",
        zip = OUTPUT_BASE + "/reads/{sample_folder_name}/02_fastqc_raw/02_{sample_r}_raw_fastqc.zip" 
        # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet",
        queue = "long" # for DKFZ LSF cluster queue 
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v5.7.0/bio/fastqc"