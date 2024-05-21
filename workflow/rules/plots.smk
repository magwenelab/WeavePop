rule depth_distribution_plots:
    input:
        rules.depth_distribution.output.distrib,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "depth_distribution.png",
        OUTDIR / "plots" / "{sample}" / "depth_global_distribution.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/depth_distribution_{sample}.log"
    script:
        "../scripts/depth_distribution_plot.R"

rule depth_by_chrom_plots:
    input:
        rules.depth_by_chrom_raw.output,
        rules.depth_by_chrom_good.output,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "depth_by_chrom.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/depth_by_chrom_{sample}.log"
    script:
        "../scripts/depth_by_chrom_plot.R"

def depth_by_regions_plot_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": OUTDIR / "mosdepth" / s["sample"]  / "depth_by_regions.tsv",
        "cnv": OUTDIR / "cnv" / s["sample"] / "cnv_calls.tsv",
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
    }
rule depth_by_regions_plot:
    input:
        unpack(depth_by_regions_plot_input),
        CHROM_NAMES,
        rules.loci.output.locitable
    output:
        OUTDIR / "plots" / "{sample}" / "depth_by_regions.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/depth_by_regions_{sample}.log"
    script:
        "../scripts/depth_by_regions_plot.R"

rule normalized_chrom_depth_plot:
    input:
        rules.dataset_metrics.output.alln
        CHROM_NAMES
    output:
        DATASET_OUTDIR / "plots" / "normalized_chromosome_depth_by_chrom.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/normalized_depth_distribution_{sample}.log"
    script:
        "../scripts/normalized_depth_distribution_plot.R"

rule dataset_summary_plot:  
    input:
        SAMPLEFILE,
        CHROM_NAMES,
        rules.dataset_metrics.output.allg,
        rules.dataset_metrics.output.allr,
        
    output:
        DATASET_OUTDIR / "plots" / "dataset_summary.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/dataset_summary.log"
    script:
        "../scripts/dataset_summary_plot.R"
