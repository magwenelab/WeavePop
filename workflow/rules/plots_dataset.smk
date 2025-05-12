# =================================================================================================
#   Per dataset | Plot dataset mapping quality and depth summary
# =================================================================================================


rule dataset_summary_plot:
    input:
        rules.quality_filter.output.metadata,
        rules.quality_filter.output.chromosomes,
        rules.join_mapping_stats.output,
    output:
        DATASET_DIR / "plots" / "dataset_summary.png",
    params:
        config["plotting"]["scale"],
    log:
        LOGS / "dataset" / "plots" / "dataset_summary.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_summary_plot.R"


# =================================================================================================
#   Per dataset | Plot normalized chromosome depth of all samples
# =================================================================================================


rule dataset_depth_by_chrom_plot:
    input:
        rules.quality_filter.output.metadata,
        rules.quality_filter.output.chromosomes,
        rules.join_cnv_chromosomes.output,
    output:
        DATASET_DIR / "plots" / "dataset_depth_by_chrom.png",
    params:
        column=config["plotting"]["metadata2color"],
        scale=config["plotting"]["scale"],
    log:
        LOGS / "dataset" / "plots" / "dataset_depth_by_chrom.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_depth_by_chrom_plot.R"
