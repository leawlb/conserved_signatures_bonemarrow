
# cutadapt for adapter (and other unwanted sequence) trimming

# using this github repository as a rough guide: 
# https://github.com/jrderuiter/snakemake-rnaseq.git

# using this this snakemake wrapper as a template but I want to be able
# to adapt this more to my needs if necessary 
# https://github.com/marcelm/cutadapt
  
# cutadapt documentation:
# https://cutadapt.readthedocs.io/en/stable/
# FOR NOW, our data has the following sequences to be cut:
# - standard illumina adapters 
# ( -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#   -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT)
# - the 3 linkers and adapter bases in read 2 remaining after umi extraction (-U 6)
# low quality ends  (-q 10)
    
rule cutadapt:
    input:
        [OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/{sample}_umiextracted_R1.fastq.gz",
         OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/{sample}_umiextracted_R2.fastq.gz"]
    output:
        fastq1 = OUTPUT_BASE + "/reads/{sample_folder_name}/04_trimmed/{sample}_trimmed_R1.fastq.gz",
        fastq2 = OUTPUT_BASE + "/reads/{sample_folder_name}/04_trimmed/{sample}_trimmed_R2.fastq.gz",
        qc =  OUTPUT_BASE + "/reads/{sample_folder_name}/04_trimmed/qc/cutadapt/{sample}.txt"
    log:
        "logs/04_cutadapt/{sample_folder_name}/{sample}.log"
    threads:
        6
    params:
        queue = "long",
        # additional options including the 
        # (adapter) read sequences to be removed are stored here:
        extra = config["cutadapt"]["extra"] 
    conda:
        "../../../envs/cutadapt1.yaml"
    resources:
        mem_mb = 5000
    shell: 
        "cutadapt "
        "--cores {threads} "
        "{params.extra} "
        "-o {output.fastq1} "
        "-p {output.fastq2} "
        "{input} > {output.qc} "
