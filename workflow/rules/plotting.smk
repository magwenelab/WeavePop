# Generate loci table
rule loci:
    input:
        refs = expand(REFDIR / "{lineage}" / "{lineage}.gff.tsv", lineage=LINEAGES),
        loci=LOCI_FILE
    output:
        locitable = DATASET_OUTDIR / "files"/ "loci_to_plot.tsv"
    log: 
        "logs/dataset/files/loci.log"
    shell:
        "xonsh workflow/scripts/loci.xsh -g {input.loci} -o {output} {input.refs} &> {log}"

# Generate coverage per chromosome plot
def coverage_plot_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "coverage": OUTDIR / "mosdepth" / s["sample"] / "good_regions_coverage.tsv",
        "structure": OUTDIR / "mosdepth" / s["sample"] / "good_structural_variants.tsv" ,
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed"),
        # "variants": DATASET_OUTDIR / "snps" / (s["group"] + "_variants.tsv")
    }
rule coverage_plot_chrom:
    input:
        unpack(coverage_plot_input),
        rules.loci.output.locitable,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "coverage.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/coverage_plot_chrom_{sample}.log"
    script:
        "../scripts/chromosome_plot.R"

# Generate coverage plots
rule coverage_stats_plot_sample:
    input:
        rules.good_coverage.output.chromosome,
        rules.raw_coverage.output.chromosome,
        CHROM_NAMES
    output:
        stats = OUTDIR / "plots" / "{sample}" / "coverage_stats.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/coverage_stats_plot_{sample}.log"
    script:
        "../scripts/coverage_stats_plots.R"

# Generate coverage stats plots
rule coverage_stats_plots_dataset:
    input:
        SAMPLEFILE,
        CHROM_NAMES,
        rules.dataset_metrics.output.allg,
        rules.dataset_metrics.output.allr,
        rules.dataset_metrics.output.allm
    output:
        DATASET_OUTDIR / "plots" / "cov_median_good.svg",
        DATASET_OUTDIR / "plots" / "cov_mean_good.svg",
        DATASET_OUTDIR / "plots" / "global.svg"
    conda:
        "../envs/r.yaml"
    params:
        config["plotting"]["metadata2color"],
        config["plotting"]["scale"]
    log:
        "logs/dataset/plots/coverage_stats_plot.log"    
    script:
        "../scripts/coverage_dataset_plots.R"

# Generate coverage distribution plots
rule cov_distribution:
    input:
        rules.samtools_stats.output.cov,
        CHROM_NAMES,
        rules.good_coverage.output.chromosome
    output:
        OUTDIR / "plots" / "{sample}" / "cov_distribution.png",
        OUTDIR / "plots" / "{sample}" / "cov_global_distribution.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/cov_distribution_{sample}.log"
    script:
        "../scripts/coverage-distribution.R"

# Generate mapq plot
def mapq_plot_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "mapq": OUTDIR / "samtools" / s["sample"] / "mapq_window.bed",
        "structure": OUTDIR / "mosdepth" / s["sample"] / "good_structural_variants.tsv" ,
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
        # "variants": DATASET_OUTDIR / "snps" / (s["group"] + "_variants.tsv")
    }
rule mapq_plot:
    input:
        unpack(mapq_plot_input),
        CHROM_NAMES,
        rules.loci.output.locitable
    output:
        OUTDIR / "plots" / "{sample}" / "mapq.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/samples/plots/mapq_plot_{sample}.log"
    script:
        "../scripts/mapq.R"

# Generate bamstats plot
# rule plot_bamstats:
#     input:
#         "analysis/{sample}/snps.bam.stats"
#     output:
#         directory("analysis/{sample}/bamstats")
#     conda:
#         "../envs/plot-bamstats.yaml"
#     log:
#         "logs/plot-bamstats/{sample}.log"
#     shell:
#         "plot-bamstats -p {output}/ {input} &> {log}"
