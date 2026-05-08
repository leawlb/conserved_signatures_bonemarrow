
# MultiQC for generating a QC summary of all kinds of QC stats

# using this snakemake wrapper: 
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/multiqc.html

#-------------------------------------------------------------------------------

# generate a list of all inputs as they do not all very nicely use the same 
# wildcards and are not all stored in the same directory

# I want to use one directory as input, but it doesn't exist until the 
# corresponding rule is exectuted. So I'm generating the directory here

# now collect all input for multiqc, including the content to be added 
# to the directory
for sn in sample_name:
    path = OUTPUT_BASE + "/alignment/" + sn + "/07_umi_tools_dedup"
    print(path)
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

# gathering all inputs together, sticking only to possible
# sample_name sample_name combinations
input_list  = []
for sn in sample_name:
    input_list = input_list + [
        expand(
            rules.fastqc_raw.output.zip,  
            sample_r = [f"{sn}_R1", f"{sn}_R2"],
            sample_folder_name = sn), 
        expand(
            rules.fastqc_prepped.output.zip,  
            sample = sn,
            sample_folder_name = sn,
            read = ["R1", "R2"]), 
        expand(
            rules.cutadapt.output.qc,  
            sample = sn,
            sample_folder_name = sn),
        expand(
            rules.star_pe_multi.output.log,
            sample = sn,
            sample_folder_name = sn),
        expand(
            rules.star_pe_multi.output.log_progress,
            sample = sn,
            sample_folder_name = sn),
        expand(
            rules.star_pe_multi.output.log_final,
            sample = sn,
            sample_folder_name = sn),
        expand( 
            OUTPUT_BASE + "/alignment/{sample_folder_name}/07_umi_tools_dedup",
            sample_folder_name = sn)
    ]

input_list = input_list + [["logs/"]]

#print(input_list)
# make the list flat, not nested
flat_list = [x for sublist in input_list for x in sublist]
print(flat_list)

# run multi qc supplying a list of files
rule test_multiqc_file:
    input:
        flat_list
    output:
        OUTPUT_BASE + "/qc/multi_qc_all/multiqc.html",
        OUTPUT_BASE + "/qc/multi_qc_all/multiqc_data.zip"
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
        queue = "long" 
    log:
        "logs/multiqc.log",
    threads: 1
    resources:
        mem_mb = 2024
    wrapper:
        "v9.4.2/bio/multiqc"
