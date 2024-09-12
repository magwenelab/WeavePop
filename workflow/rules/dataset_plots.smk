
# =================================================================================================
#   Per dataset | Plot dataset mapping quality and depth summary
# =================================================================================================

rule dataset_summary_plot:  
    input:
        INT_DATASET_DIR / "metadata.csv",
        rules.copy_config.output.c,
        rules.join_depth_by_chrom_good.output,
        rules.join_depth_by_chrom_raw.output,
        rules.join_mapping_stats.output
    output:
        DATASET_DIR / "plots" / "dataset_summary.png"
    conda:
        "../envs/r.yaml"
    params:
        config["plotting"]["scale"]
    log:
        "logs/dataset/plots/dataset_summary.log"
    script:
        "../scripts/dataset_summary_plot.R"

# =================================================================================================
#   Per dataset | Plot normalized chromosome depth of all samples
# =================================================================================================

rule dataset_depth_by_chrom_plot:
    input:
        INT_DATASET_DIR / "metadata.csv",
        rules.copy_config.output.c,
        rules.join_depth_by_chrom_good.output
    output:
        DATASET_DIR / "plots" / "dataset_depth_by_chrom.png"
    conda:
        "../envs/r.yaml"
    params:
        column = config["plotting"]["metadata2color"],
        scale = config["plotting"]["scale"]
    log:
        "logs/dataset/plots/dataset_depth_by_chrom.log"
    script:
        "../scripts/dataset_depth_by_chrom_plot.R"