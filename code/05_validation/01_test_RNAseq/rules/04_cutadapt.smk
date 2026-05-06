
# cutadapt for adapter (and other unwanted sequence) trimming

# using this github repository roughly as a guide: 
# https://github.com/jrderuiter/snakemake-rnaseq.git

# using this this snakemake wrapper as a template but I want to be able
# to adapt this more to my needs if necessary 
# https://github.com/marcelm/cutadapt
  
# cutadapt documentation:
# https://cutadapt.readthedocs.io/en/stable/
    
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
        queue = "medium",
        # additional options including the 
        # (adapter) read sequences to be removed are stored here:
        extra = config["cutadapt"]["extra"] 
    conda:
        "../../../envs/cutadapt1.yaml"
    resources:
        mem_mb = 2000
    shell: 
        "cutadapt "
        "--cores {threads} "
        "{params.extra} "
        "-o {output.fastq1} "
        "-p {output.fastq2} "
        "{input} > {output.qc} "

### explanation of option choices (see also config for extra parameters) ### 

### UMI linker: -U 6
# "Trim 8 nt UMIs + 3 nt UMI linker + 3 nt from the SMART UMI Adapter from Read2 prior to mapping" 
# from the SMART-Seq® Total RNA Library Prep with ZapR Depletion (withUMIs) User Manual
# the UMI has already been removed by umi_tools, so 6 nt remain

### read length is longer than insert size: trimming the 3' end of reads containing the rest of the adapter sequence
### -a and -A
# https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314
# Truseq adapters: https://dnatech.ucdavis.edu/sites/g/files/dgvnsk15531/files/inline-files/illumina-adapter-sequences_1000000002694-00.pdf
# mentioned also in cutadapt: https://cutadapt.readthedocs.io/en/stable/guide.html#truseq
# cutadapt will also remove the sequence following the adapters

### low quality ends: -q N
# trims low quality below quality N
# because the quality according to fastQC remains high, no quality trimming for now

### low quality ends including the detection of "high-quality Gs" from what was actually dark circles on the machine
### --nextseq-trim N

### poly-A
# no strong indication in fastQC of a lot of TTT or AAA so far, so also no poly-A trimming for now
