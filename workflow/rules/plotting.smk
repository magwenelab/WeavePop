# Generate loci table
rule loci:
    input:
        refs = expand(rules.gff2tsv.output, lineage=LINEAGES),
        loci=LOCI_FILE
    output:
        locitable = DATASET_OUTDIR / "files"/ "loci_to_plot.tsv"
    log: 
        "logs/references/loci.log"
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
        OUTDIR / "plots" / "{sample}" / "coverage.svg"
    conda:
        "../envs/r.yaml"
    log:
        "logs/coverage/coverage_plot_chrom_{sample}.log"
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
        "logs/coverage/coverage_stats_plot_{sample}.log"
    script:
        "../scripts/coverage_stats_plots.R"

# Generate coverage stats plots
rule coverage_stats_plots_dataset:
    input:
        SAMPLEFILE,
        rules.dataset_metrics.output.allg,
        rules.dataset_metrics.output.allr,
        CHROM_NAMES
    output:
        DATASET_OUTDIR / "plots" / "cov_median_good.png",
        DATASET_OUTDIR / "plots" / "cov_mean_good.png",
        DATASET_OUTDIR / "plots" / "cov_global.png"
    conda:
        "../envs/r.yaml"
    params:
        config["plotting"]["metadata2color"]
    log:
        "logs/coverage/stats_plot.log"    
    script:
        "../scripts/coverage_dataset_plots.R"

# Generate mapq distribution plot
rule mapq_distribution:
    input:
        rules.samtools_stats.output.mapq,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "mapq_distribution.png"
    conda:
        "../envs/r.yaml"   
    log:
        "logs/mapq/mapq_distribution_{sample}.log"
    script:
        "../scripts/mapq-distribution.R"

# Generate coverage distribution plots
rule cov_distribution:
    input:
        rules.samtools_stats.output.cov,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "cov_distribution.png",
        OUTDIR / "plots" / "{sample}" / "cov_global_distribution.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/coverage/cov_distribution_{sample}.log"
    script:
        "../scripts/coverage-distribution.R"

# Generate mapped reads plot
rule mapped_plot:
    input:
        rules.dataset_metrics.output.allm,
        SAMPLEFILE
    output:
        DATASET_OUTDIR / "plots" / "mapped_reads.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/stats/mapped_plot.log"
    script:
        "../scripts/mapped_reads.R"

# Generate mapq plot
rule mapq_plot:
    input:
        rules.mapq.output.winbed,
        CHROM_NAMES,
        rules.loci.output.locitable
    output:
        OUTDIR / "plots" / "{sample}" / "mapq.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/mapq/mapq_plot_{sample}.log"
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
