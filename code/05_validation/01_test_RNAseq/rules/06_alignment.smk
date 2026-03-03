
# STAR for alignment of pre-processed fastq files

# using star snakemake wrapper 
# https://snakemake-wrappers.readthedocs.io/en/v3.5.3/wrappers/star/align.html

# STAR documentation
# https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf

# first, generate STAR genome index
# I will download the current mus musculus genome from ensembl, release 115
# even though it is not the same used for scRNAseq alignment
rule download_genome:
    output:
        fasta = OUTPUT_BASE + "/genome/fasta/Mus_musculus.GRCm39.dna.primary_assembly.fa",
        gtf = OUTPUT_BASE + "/genome/gtf/Mus_musculus.GRCm39.115.gtf"
    params:
        url_fasta = "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz",
        url_gtf = "https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz",
        genome_directory = OUTPUT_BASE + "/genome"
    shell:
        r"""
        #cd {params.genome_directory};
        #pwd;

        echo {params.url_fasta}
        echo {params.url_gtf}

        # these are absolute paths
        mkdir -p fasta;
        wget -O {output.fasta}.gz {params.url_fasta};
        gunzip -k -f {output.fasta}.gz;

        mkdir -p gtf;
        wget -O {output.gtf}.gz {params.url_gtf};
        gunzip -k -f {output.gtf}.gz;

        """

rule generate_star_index:
    input:
        fasta = rules.download_genome.output.fasta,
        gtf = rules.download_genome.output.gtf
    output:
        index_dir = directory(OUTPUT_BASE + "/genome/star")
    threads:
        12
    params:
        STAR_PATH = config["star_path"],
        queue = "medium",
        mem_mb = 40000
    shell:
        r"""
        mkdir -p {output.index_dir};

        "{params.STAR_PATH}" \ 
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 49 
        """

# invoke star using the absolute path to work around module problems
# sjdbOverhang should be number of reads -1, so 50bp -1 = 49