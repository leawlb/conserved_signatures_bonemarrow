
# umi extraction to allow for adapter trimming alignment without UMIs

# umi_tools documentation:
# https://umi-tools.readthedocs.io/en/latest/QUICK_START.html
# "UMIs are removed and appended to the read name"
# "Any other barcode, for example a library barcode, is left on the read"

# nf-core/rnaseq also does umi extraction first
# https://github.com/shahcompbio/bulk-illumina-rnaseq.git

rule umi_tools_extract:
    input:
        # raw files (symlinks)
        fastq1 = OUTPUT_BASE + "/reads/{sample_folder_name}/01_raw/01_{sample}_R1_raw.fastq.gz",
        fastq2 = OUTPUT_BASE + "/reads/{sample_folder_name}/01_raw/01_{sample}_R2_raw.fastq.gz"
    output: 
        fastq1 = OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/03_{sample}_umiextracted_R1.fastq.gz",
        fastq2 = OUTPUT_BASE + "/reads/{sample_folder_name}/03_umis_extracted/03_{sample}_umiextracted_R2.fastq.gz"
    log:
        "logs/{sample_folder_name}/03_umi_tools_extract/03_{sample}.log"
    threads:
        1
    params:
        queue = "long",
        extra = config["umi_tools_extract"]["extra"] 
    conda:
        "../../../envs/umi_tools5.yaml"
    resources:
        mem_mb = 2000
    shell:
        "umi_tools extract "
        "--extract-method=string "
        "{params.extra} "
        "--read2-in={input.fastq1} "
        "--read2-out={output.fastq1} "
        "-L {log} "
        "-I {input.fastq2} "
        "-S {output.fastq2} "

### explanation of option choices (see also config for extra parameters) ### 
# "Trim 8 nt UMIs + 3 nt UMI linker + 3 nt from the SMART UMI Adapter from Read2 prior to mapping" 
# according to SMART-Seq® Total RNA Library Prep with ZapR Depletion (with UMIs) User Manual

# --bc-pattern:
# NNNNNNNN will be added to the sequence name (that is the UMI) from the start of the read
# XXXXXX will be kept in the read, then I will trim it downstream using cutadapt

# UMI is on read 2 only; it is recommended to swap the input and output files from R1 to R2
# - https://github.com/CGATOxford/UMI-tools/issues/331
# - https://github.com/CGATOxford/UMI-tools/issues/522
# there is a --read2-only option but it's not available for this version of umi_tools
# otherwise, umi_tools also expects both -bc-pattern and bc-pattern2
