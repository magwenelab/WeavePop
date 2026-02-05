# =================================================================================================
#   Per dataset | Plot dataset mapping quality and depth summary
# =================================================================================================


rule dataset_summary_plot:
    input:
        metadata=rules.quality_filter.output.metadata,
        chroms=rules.quality_filter.output.chromosomes,
        stats=rules.join_mapping_stats.output,
    output:
        plot=DATASET_DIR / "plots" / "mapping_summary.png",
    params:
        scale=config["plotting"]["scale"],
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
        metadata=rules.quality_filter.output.metadata,
        chroms=rules.quality_filter.output.chromosomes,
        cnv=rules.join_cnv_chromosomes.output,
    output:
        plot=DATASET_DIR / "plots" / "dataset_depth_by_chrom.png",
    params:
        column=COLOR_BY,
        scale=config["plotting"]["scale"],
    log:
        LOGS / "dataset" / "plots" / "dataset_depth_by_chrom.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_depth_by_chrom_plot.R"


# =================================================================================================
#   Per ref_genome | Plot classification of variants
# =================================================================================================

rule refs_variant_classification_plots:
    input:
        chromosomes=DATASET_DIR / "chromosomes.csv",
        classif=INT_DATASET_DIR / "snpeff" / "{ref_genome}_variant_classification.tsv",
        variants=INT_DATASET_DIR / "snpeff" / "{ref_genome}_variants.tsv",
        presence=INT_DATASET_DIR / "snpeff" / "{ref_genome}_presence.tsv",
        loci=INT_REFS_DIR / "{ref_genome}" / "loci_to_plot.tsv",
    output:
        plot=DATASET_DIR / "plots" / "{ref_genome}_variant_summary.png",
        plot_density=REFS_DIR / "{ref_genome}" / "{ref_genome}_variants_by_windows.png",
    params:
        window_size=config["depth_quality"]["mosdepth"]["window"],
    log:
        LOGS / "dataset" / "plots" / "refs_variant_classification_plots_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/refs_variant_classification_plots.R"