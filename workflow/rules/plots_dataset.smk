# =================================================================================================
#   Per dataset | Summary of mapping stats and variant calling 
# =================================================================================================


rule dataset_summary_plot:
    input:
        metadata=rules.quality_filter.output.metadata,
        chroms=rules.quality_filter.output.chromosomes,
        stats=rules.join_mapping_stats.output,
    output:
        plot=DATASET_DIR / "plots" / "mapping_summary.png",
    log:
        LOGS / "dataset" / "plots" / "dataset_summary.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_summary_plot.R"


# =================================================================================================
#   Per dataset | Summary plots of depth, cnv and variants
# =================================================================================================


rule dataset_depth_plot:
    input:
        metadata=rules.quality_filter.output.metadata,
        chroms=rules.quality_filter.output.chromosomes,
        cnv=rules.join_cnv_chromosomes.output,
    output:
        plot=DATASET_DIR / "plots" / "depth_summary.png",
    params:
        column=COLOR_BY,
    log:
        LOGS / "dataset" / "plots" / "dataset_depth_plot.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_depth_plot.R"


rule dataset_variants_plot:
    input:
        metadata=rules.quality_filter.output.metadata,
        classif=DATASET_DIR / "snpeff" / "variant_classification.tsv",
        variants=DATASET_DIR / "snpeff" / "variants.tsv",
        presence=DATASET_DIR / "snpeff" / "presence.tsv",
    output:
        plot=DATASET_DIR / "plots" / "variant_summary.png",
    log:
        LOGS
        / "dataset"
        / "plots"
        / "dataset_variants_plot.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_variants_plot.R"


rule dataset_cnv_plot:
    input:
        metadata=rules.quality_filter.output.metadata,
        chromosomes=DATASET_DIR / "chromosomes.csv",
        cnv=DATASET_DIR / "cnv" / "cnv_chromosomes.tsv",
    output:
        plot=DATASET_DIR / "plots" / "cnv_summary.png",
        plot2=DATASET_DIR / "plots" / "depth_uniformity_summary.png",
    log:
        LOGS
        / "dataset"
        / "plots"
        / "dataset_cnv_plot.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/dataset_cnv_plot.R"

# =================================================================================================
#   Per ref_genome | Plot of classification of variants
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
        LOGS
        / "dataset"
        / "plots"
        / "refs_variant_classification_plots_{ref_genome}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/refs_variant_classification_plots.R"