
# generate useful symbolic links for each fastq file

#-------------------------------------------------------------------------------

# this is specific for DKFZ-type output and folder structure since I'm directly
# using the DKFZ metadata tsv

# create a sample-specific directory for each sample first
# then the rule can simply add both symlinks (both reads) into each directory
for sn in sample_name:
    dir_path = OUTPUT_BASE + "/reads/" + sn + "/01_raw"  
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)

# for each sample, there are two fastq files, one R1 and one R2
# this rule generates one symlink per R1/R2 file and gives interpretable names
rule link_files:
    output: 
        OUTPUT_BASE + "/reads/{sample_folder_name}/01_raw/01_{sample_r}_raw.fastq.gz"
    params:
        queue = "short" # for DKFZ LSF cluster queue 
    threads: 1
    resources:
        mem_mb = 500
    run:
        print(output[0])
        print(wildcards.sample_r)

        # for each sample (with R1/R2 info), find the associated fastq file name 
        # using the metadata table
        fastq_file_path =  metadata_bulk.loc[metadata_bulk["SAMPLE_NAME_R"] == wildcards.sample_r, "FastQ Path"].iloc[0]

        print(fastq_file_path)

        # generate symbolic link using the sample name to the original
        # fastq file location
        if os.path.exists(output[0]):
            os.remove(output[0])
        os.symlink(fastq_file_path, output[0])
