# =================================================================================================
#   Join all lineages | Create table with loci to add to plots
# =================================================================================================


rule loci:
    input:
        gff=REFS_DIR / "{lineage}" / "{lineage}.gff.tsv",
        loci=LOCI_FILE,
    output:
        locitable=INT_REFS_DIR / "{lineage}" / "loci_to_plot.tsv",
    log:
        LOGS / "references" / "plots" / "loci_{lineage}.log",
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/loci.py"


# =================================================================================================
#   Per sample | Plot depth distribution and depth by chromosome
# =================================================================================================


rule depth_distribution_plots:
    input:
        distrib=INT_SAMPLES_DIR
        / "depth_quality"
        / "{sample}"
        / "depth_distribution.tsv",
        chroms=CHROM_NAMES_FILE,
        metadata=METADATA_ORIGINAL_FILE,
    output:
        chrom=SAMPLES_DIR / "plots" / "{sample}" / "depth_chrom_distribution.png",
        globl=SAMPLES_DIR / "plots" / "{sample}" / "depth_global_distribution.png",
    conda:
        "../envs/r.yaml"
    log:
        LOGS / "samples" / "plots" / "depth_distribution_{sample}.log",
    script:
        "../scripts/depth_distribution_plots.R"


# =================================================================================================
#   Per sample | Plot depth and mapq by windows
# =================================================================================================


rule depth_boxplot:
    input:
        unpack(depth_boxplot_input),
        metadata=METADATA_ORIGINAL_FILE,
    output:
        plot=SAMPLES_DIR / "plots" / "{sample}" / "depth_boxplot.png",
    conda:
        "../envs/r.yaml"
    log:
        LOGS / "samples" / "plots" / "depth_boxplot_{sample}.log",
    script:
        "../scripts/depth_boxplots.R"


rule depth_by_windows_plots:
    input:
        unpack(depth_by_windows_plots_input),
        metadata=METADATA_ORIGINAL_FILE,
    output:
        plot=SAMPLES_DIR / "plots" / "{sample}" / "depth_by_windows.png",
    params:
        find_repeats=config["repeats"]["activate"],
    conda:
        "../envs/r.yaml"
    log:
        LOGS / "samples" / "plots" / "depth_by_windows_{sample}.log",
    script:
        "../scripts/depth_by_windows_plots.R"


rule depth_vs_cnvs_plots:
    input:
        unpack(depth_vs_cnvs_plots_input),
        metadata=METADATA_ORIGINAL_FILE,
    output:
        SAMPLES_DIR / "plots" / "{sample}" / "depth_vs_cnvs.png",
    conda:
        "../envs/r.yaml"
    log:
        LOGS / "samples" / "plots" / "depth_vs_cnvs_{sample}.log",
    script:
        "../scripts/depth_vs_cnvs_plots.R"


rule mapq_plots:
    input:
        unpack(mapq_plot_input),
        metadata=METADATA_ORIGINAL_FILE,
    output:
        SAMPLES_DIR / "plots" / "{sample}" / "mapq.png",
    params:
        find_repeats=config["repeats"]["activate"],
    conda:
        "../envs/r.yaml"
    log:
        LOGS / "samples" / "plots" / "mapq_plots_{sample}.log",
    script:
        "../scripts/mapq_plots.R"


        
rule variant_classification_plots:
    input:
        unpack(variant_classification_plots_input)
    output:
        barplot=SAMPLES_DIR / "plots" / "{sample}" / "variants_barplot.png",
        status=SAMPLES_DIR / "plots" / "{sample}" / "variants_by_windows_status.png",
        impact=SAMPLES_DIR / "plots" / "{sample}" / "variants_by_windows_impact.png",
    params:
        window_size=config["depth_quality"]["mosdepth"]["window"],
    conda:
        "../envs/r.yaml"
    log:
        LOGS / "samples" / "plots" / "variant_classification_plots_{sample}.log",
    script:
        "../scripts/variant_classification_plots.R"
