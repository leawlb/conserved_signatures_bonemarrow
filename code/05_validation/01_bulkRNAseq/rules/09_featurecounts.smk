
# FeatureCounts to get reads per gene

# loosely based on this snakemake wrapper:
# https://snakemake-wrappers.readthedocs.io/en/0.72.0/wrappers/subread/featurecounts.html
# and this repository:
# https://github.com/jrderuiter/snakemake-rnaseq/blob/develop/wrappers/subread/feature_counts/wrapper.py

# FeatureCounts documentation
# https://subread.sourceforge.net/featureCounts.html

#-------------------------------------------------------------------------------

# get all input bam files, which are the depuplicated, sorted, indexed ones
input_list = []
for sn in sample_name:
    input_list = input_list + [
        expand(
            rules.generate_sorted_dedup_bam.output.sorted_dedup_bam,  
            sample = sn,
            sample_folder_name = sn)]

flat_list = [x for sublist in input_list for x in sublist]
#print(flat_list)

rule feature_counts:
    input:
        sample_list = flat_list, # list of sam or bam files
        annotation = rules.download_genome.output.gtf
    output:
        featurecounts = OUTPUT_BASE + "/counts/09_feature_counts/09_all_counts.featureCounts",
        summary = OUTPUT_BASE + "/counts/09_feature_counts/09_all_counts.featureCounts.summary",
        jcounts = OUTPUT_BASE + "/counts/09_feature_counts/09_all_counts.featureCounts.jcounts"
    threads:
        2
    params:
        extra = config["feature_counts"]["extra"],
        queue = "long"
    resources:
        mem_mb = 2000
    log:
        "logs/counts/09_featurecounts/09_all_counts.log"
    conda:
        "../../../envs/feature_counts.yaml"
    shell:
        "featureCounts "
        " -a {input.annotation} "
        " {params.extra} "
        " -o {output.featurecounts} "
        " -T {threads} "
        " {input.sample_list} "

# option choices

# -p: paired end input
# --countReadPairs_ count reads instead of fragments, because PE
# -t exon: use exon type features from annotation
# -g gene_id: attribute used to group features = count at gene level
# -J: from wrapper, related to splicing junctions
# -O: reads overlapping multiple features are counted in all features
# --fracOverlap 0.2: at least 20% of a read must overlap with a feature to count