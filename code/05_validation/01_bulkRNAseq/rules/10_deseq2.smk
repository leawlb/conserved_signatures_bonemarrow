
# read counts into DESeq2 object, check some qc, perform analysis to 
# compare per target between conditions to figure out which populations
# were isolated

#-------------------------------------------------------------------------------

print(config["base"] + config["metadata_paths"]["ensembl_mus"])
print(rules.feature_counts.output.featurecounts)
print(config["base"] + config["bulk_rnaseq_data_paths"])

# read counts and metadata into a DESeq2 object (dds)
rule read_dds:
    input:
        ensembl_mus_path = config["base"] + config["metadata_paths"]["ensembl_mus"],
        counts_path = rules.feature_counts.output.featurecounts,
        metadata_path = config["table_bulk"]
    params:
        queue="short"
    resources:
        mem_mb=2000,
    threads: 
        2
    output:
        output_path_dds = OUTPUT_BASE + "/counts/10_dds/10_dds.Rds"
    conda:
        "../../../envs/run_deseq2_val.yaml"
    script:
        "../scripts/10_read_dds.R"

print( OUTPUT_BASE + "/counts/11_deseq2/11_vsd_{ab_target}.Rds")
print(OUTPUT_BASE + "/counts/11_deseq2/11_res_df_{ab_target}.Rds")

# qc and analysis using DESeq2
# also separate dds object by antibody target 
rule deseq2_qc:
    input:
        input_path_dds = rules.read_dds.output.output_path_dds
    params:
        queue="short"
    resources:
        mem_mb=2000,
    threads: 
        2
    output:
        output_path_vsd = OUTPUT_BASE + "/counts/11_deseq2/11_vsd_{ab_target}.Rds",
        output_path_results = OUTPUT_BASE + "/counts/11_deseq2/11_res_df_{ab_target}.Rds"
    conda:
        "../../../envs/run_deseq2_val.yaml"
    script:
        "../scripts/11_deseq_qc.R"

# run report on DEseq2 quality
rule check_deseq2_qc_report:
    input:
        dds_path = rules.read_dds.output.output_path_dds,
        vsd_path = rules.deseq2_qc.output.output_path_vsd,
        results_df_path = rules.deseq2_qc.output.output_path_results
    params:
        queue="medium"
    resources:
        mem_mb=4000,
    threads: 
        1
    output:
        OUTPUT_BASE + "/counts/11_deseq2/11_check_deseq2_qc_{ab_target}.html"
    conda:
        "../../../envs/run_deseq2_val.yaml"
    script:
        "../check_deseq2_qc.Rmd"

# run figure_ideas for each antibody target, generating multiple plots
rule figure_ideas_neo1:
    input:
        vsd_neo1 = expand(rules.deseq2_qc.output.output_path_vsd, ab_target = "neo1"),
        results_df_neo1 = expand(rules.deseq2_qc.output.output_path_results, ab_target = "neo1"),
        sce_str = config["base"] + config["scRNAseq_data_paths"]["main"] + "/sce_objects/02_sce_anno/10_anns/sce_str.rds",
        sign_genes = config["base"] + config["scRNAseq_data_paths"]["main"] + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own/01_gens/geneset_list_str"
    params:
        queue="medium",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    resources:
        mem_mb=40000,
    threads: 
        1
    output:
        OUTPUT_BASE + "/counts/11_deseq2/11_figure_ideas_neo1.html"
    conda:
        "../../../envs/bulk_figures5.yaml"
    script:
        "../figure_ideas_neo1.Rmd"

rule figure_ideas_trkb:
    input:
        vsd_trkb = expand(rules.deseq2_qc.output.output_path_vsd, ab_target = "trkb"),
        results_df_trkb = expand(rules.deseq2_qc.output.output_path_results, ab_target = "trkb"),
        sce_str = config["base"] + config["scRNAseq_data_paths"]["main"] + "/sce_objects/02_sce_anno/10_anns/sce_str.rds",
        sign_genes = config["base"] + config["scRNAseq_data_paths"]["main"] + "/sce_objects/03_sce_analysis/04_signatures/01_reclustering_own/01_gens/geneset_list_str"
    params:
        queue="medium",
        colors_path = COLORS,
        colors = "../../source/colors.R",
        plotting = "../../source/plotting.R"
    resources:
        mem_mb=40000,
    threads: 
        1
    output:
        OUTPUT_BASE + "/counts/11_deseq2/11_figure_ideas_trkb.html"
    conda:
        "../../../envs/bulk_figures5.yaml"
    script:  
        "../figure_ideas_trkb.Rmd"