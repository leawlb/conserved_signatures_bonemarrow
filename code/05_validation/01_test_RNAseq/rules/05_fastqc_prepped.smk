
# FASTQC for generating a QC html for each fastq file

# check quality on individual fastq files after some pre-processing
# using this snakemake wrapper:
# https://github.com/s-andrews/FastQC

# nf-core/rnaseq also does FastQC step at this point
# https://github.com/shahcompbio/bulk-illumina-rnaseq.git

# because the prior rule runs per sample and has one output item for each read 
# (fastq1 and fastq2)
# but I want to perform FastQC one each read separately

def cutadapt_output_fastq(wc):
    #print(wc.read)
    if wc.read == "R1":
        return rules.cutadapt.output.fastq1
    else:
        return rules.cutadapt.output.fastq2

# depending on the wildcard used for each fastqc_prepped instance,
# either fastq1 or fastq2 cutadapt output item will be used as input 

rule fastqc_prepped:
    input:
        cutadapt_output_fastq
    output:
        html = OUTPUT_BASE + "/reads/{sample_folder_name}/05_fastqc_prepped/{sample}_{read}_prepped.html",
        zip = OUTPUT_BASE + "/reads/{sample_folder_name}/05_fastqc_prepped/{sample}_{read}_prepped_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet ",
        queue = "short"  
    log:
        "logs/05_fastqc_prepped/{sample_folder_name}/{sample}_{read}.log"
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v5.7.0/bio/fastqc"

