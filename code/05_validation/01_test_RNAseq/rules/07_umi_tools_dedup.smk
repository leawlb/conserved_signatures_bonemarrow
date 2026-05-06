
# umi deduplication after alignment to mark duplicates due to PCR amplification

# umi_tools documentation:
# https://umi-tools.readthedocs.io/en/latest/QUICK_START.html
# "UMIs are removed and appended to the read name"
# "Any other barcode, for example a library barcode, is left on the read"

# nf-core/rnaseq also does umi dedup after alignment
# https://nf-co.re/rnaseq/3.14.0/

#-------------------------------------------------------------------------------

# umi_tools dedup uses .bam, follow their example for generating suitable files 
rule generate_bam:
    input:
        # I use "aligned" for consistency, umi_tools uses "mapped"
        aligned_sam = rules.star_pe_multi.output.aln
    output:
        aligned_sorted_bam = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_sorted_aligned.bam",
        aligned_sorted_bam_index = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_sorted_aligned.bam.bai",
    log:
        "logs/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_generate_bam.log"   
    threads:
        1
    params:
        queue = "long",
        output_dir = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup"
    envmodules:
        "SAMtools/1.20-GCC-14.1.0"
    resources:
        mem_mb = 2000
    shell:
        "set -euo pipefail; "
        "mkdir -p {params.output_dir}; "
        "samtools view -bS {input.aligned_sam} > {params.output_dir}/aligned.bam; "
        "samtools sort {params.output_dir}/aligned.bam -o {output.aligned_sorted_bam}; "
        "samtools index {output.aligned_sorted_bam} {output.aligned_sorted_bam_index}; "

# get statistics on the dedupped sorted bam file for QC BEFORE dedup
rule samtools_stats_before:
    input:
        aligned_sorted_bam = rules.generate_bam.output.aligned_sorted_bam
    output:
        samtools_stats_before = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_before_stats.txt"
    log:
        "logs/{sample_folder_name}/07_umi_tools_dedup/07_before_{sample}_samtools_stats.log"   
    threads:
        1
    params:
        queue = "medium",
    resources:
        mem_mb = 500
    envmodules:
        "SAMtools/1.20-GCC-14.1.0"
    shell:
        "set -euo pipefail; "
        "samtools stats {input.aligned_sorted_bam} > {output.samtools_stats_before} "

# use umi_tools dedup function to deduplicate the bam file
rule umi_tools_dedup:
    input:
        aligned_sorted_bam = rules.generate_bam.output.aligned_sorted_bam
    output: 
        dedup_bam = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_dedup.bam",
    log:
        "logs/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_dedup.log"   
    threads:
        1
    params:
        queue = "long",
        output_dir = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup"
    conda:
        "../../../envs/umi_tools5.yaml"
    resources:
        mem_mb = 8000
    shell:
        "umi_tools dedup -I {input.aligned_sorted_bam} --log={log} --paired --output-stats={params.output_dir}/deduplicated --stdout {output.dedup_bam} "

# from dedup output, generate sorted and indexed .bam for downstream processing
# this rule also generates some QC outputs in the directory but that are not
# explicitly stated in output
# this takes ~30 - 40 minutes
rule generate_sorted_dedup_bam:
    input:
        dedup_bam = rules.umi_tools_dedup.output.dedup_bam
    output:
        sorted_dedup_bam = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_sorted_dedup.bam",
        sorted_dedup_bam_index = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_sorted_dedup.bam.bai"
    log:
        "logs/07_umi_tools_dedup/{sample_folder_name}/{sample}_generate_bam.log"   
    threads:
        1
    params:
        queue = "long",
        output_dir = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup"
    envmodules:
        "SAMtools/1.20-GCC-14.1.0"
    resources:
        mem_mb = 2000
    shell:
        "set -euo pipefail; "
        "samtools sort {input.dedup_bam} -o {output.sorted_dedup_bam}; "
        "samtools index {output.sorted_dedup_bam} {output.sorted_dedup_bam_index}; "

# get statistics on the dedupped sorted bam file for QC AFTER dedup
rule samtools_stats_dedup: 
    input:
        sorted_dedup_bam = rules.generate_sorted_dedup_bam.output.sorted_dedup_bam
    output:
        samtools_stats_after = OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup/07_{sample}_dedup_stats.txt"
    log:
        "logs/{sample_folder_name}/07_umi_tools_dedup/07_dedup_{sample}_samtools_stats.log"   
    threads:
        1
    params:
        queue = "medium",
    resources:
        mem_mb = 500
    envmodules:
        "SAMtools/1.20-GCC-14.1.0"
    shell:
        "set -euo pipefail; "
        "samtools stats {input.sorted_dedup_bam} > {output.samtools_stats_after} "
