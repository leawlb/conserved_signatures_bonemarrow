
# specific for DKFZ-type output and folder structure
# using the DKFZ metadata_tsv
# check regex when implementing on new files

rule link_files:
    output: 
        OUTPUT_BASE + "/reads/{sample_folder_name}/{sample_r}.fastq.gz"
    run:
        #print(output[0])
        #print(wildcards.sample_r)

        # for each sample (with R1/R2 info), find the associated fastq file name 
        # using the metadata table
        fastq_file_name = metadata_bulk.loc[metadata_bulk["SAMPLE_NAME_R"] == wildcards.sample_r, "FASTQ_FILE"].iloc[0]
        fastq_folder_name = re.sub(r'_R[12]\.fastq\.gz$', '', fastq_file_name)
        fastq_file_path = FASTQ_INPUT_DIR + "/" + fastq_folder_name + "/fastq/" + fastq_file_name

        #print(fastq_file_name)
        #print(fastq_folder_name)
        print(fastq_file_path)

        # generate symbolic link using the sample name to the original
        # fastq file location
        if os.path.exists(output[0]):
            os.remove(output[0])
        os.symlink(fastq_file_path, output[0])


