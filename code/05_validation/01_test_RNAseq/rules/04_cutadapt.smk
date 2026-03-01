
# cutadapt for adapter (and other unwanted sequence) trimming

# using this github repository as a rough guide: 
# https://github.com/jrderuiter/snakemake-rnaseq.git

# using this this snakemake wrapper as a template but I want to be able
# to adapt this to more to my needs if necessary 
# https://github.com/marcelm/cutadapt
  
# cutadapt documentation:
# https://cutadapt.readthedocs.io/en/stable/
# our data has the following sequences to be cut:
# - TODO: FILL IN
    
rule cutadapt:
    input:
        [OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/{sample}_R1.fastq.gz",
         OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/{sample}_R2.fastq.gz"]
    output:
        fastq1 = OUTPUT_BASE + "/reads/{sample_folder_name}/04_trimmed/{sample}_R1.fastq.gz"
        fastq2 = OUTPUT_BASE + "/reads/{sample_folder_name}/04_trimmed/{sample}_R2.fastq.gz"
        qc =  OUTPUT_BASE + "/reads/{sample_folder_name}/04_trimmed/qc/cutadapt/{sample}.txt"
    log:
        "logs/04_cutadapt/{sample_folder_name}/{sample}.log"
    threads:
        6
    params:
        queue = "long",
        # additional options including the 
        # (adapter) read sequences to be removed are stored here:
        config["cutadapt"]["extra"] 
    conda:
        "../../envs/cutadapt.yaml"
    shell: 
        "cutadapt "
        "--cores {threads} "
        "{params.extra} "
        "-o {output.fastq1} " # trimmed R1 
        "-p {output.fastq2} " # trimmed R2
        "{input} > {output.qc} {log}"
