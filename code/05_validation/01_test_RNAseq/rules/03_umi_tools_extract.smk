
# umi extraction to allow for adapter trimming alignment without UMIs

# umi_tools documentation:
# https://umi-tools.readthedocs.io/en/latest/QUICK_START.html
# "UMIs are removed and appended to the read name"
# "Any other barcode, for example a library barcode, is left on the read"

# nf-core/rnaseq also does umi extraction first
# https://github.com/shahcompbio/bulk-illumina-rnaseq.git

rule: umi_tools_extract:
    input:
        # raw files (symlinks)
        fastq1 = OUTPUT_BASE + "/reads/{sample_folder_name}/{sample}_R1.fastq.gz"
        fastq2 = OUTPUT_BASE + "/reads/{sample_folder_name}/{sample}_R2.fastq.gz"
    output: 
        fastq1 = OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/{sample}_R1.fastq.gz"
        fastq2 = OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/{sample}_R2.fastq.gz"
    log:
        "logs/03_umi_tools_extract/{sample_folder_name}/{sample}.log"
    threads:
        1
    params:
        queue = "short",
        # additional options including the UMI barcode info is stored here:
        config["umi_tools_extract"]["extra"] 
    conda:
        "../../envs/umi_tools.yaml"
    shell:
        "umi_tools extract "
        "--extract-method=string "
        "{params.extra} "
        "--read2-in={input.fastq2} "
        "--read2-out={output.fastq2} "
        "-L {log} "
        "{input.fastq1} > {output.fastq1} "
