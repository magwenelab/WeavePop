configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))

rule all:
    input:
        expand("analysis/{sample}/coverage.svg",sample=samples),
        "results/cov_median_good.svg",
        expand("analysis/{sample}/mapq_distribution.svg",sample=samples),
        expand("analysis/{sample}/cov_distribution.svg",sample=samples),
        # expand("analysis/{sample}/bamstats", sample=samples),
        "results/mapped_reads.svg",
        expand("analysis/{sample}/mapq.svg", sample=samples)

rule coverage_plot:
    input:
        "analysis/{sample}/coverage.regions.bed.gz",
        "analysis/{sample}/coverage_good.regions.bed.gz",
        "files/chromosome_names.csv",
        "files/loci_to_plot.tsv"
    output:
        "analysis/{sample}/coverage.svg",
        "analysis/{sample}/coverage_stats.svg",
        "analysis/{sample}/coverage_stats_good.csv",
        "analysis/{sample}/coverage_stats_raw.csv"
    conda:
        "envs/r.yaml"
    log:
        "logs/coverage/{sample}.log"
    script:
        "scripts/coverage.R"

rule cat_stats:
    input:
        r = expand("analysis/{sample}/coverage_stats_raw.csv",sample=samples),
        g = expand("analysis/{sample}/coverage_stats_good.csv",sample=samples),
    output:
        allr = "results/coverage_raw.csv",
        allg = "results/coverage_good.csv",
    log:
        "logs/coverage/cat_stats.log"
    shell:
        "cat {input.r} | grep -v Sample > {output.allr} "
        "&& "
        "cat {input.g} | grep -v Sample > {output.allg} "
        "2> {log}"

rule coverage_stats_plots:
    input:
        config["sample_file"],
        "results/coverage_good.csv",
        "results/coverage_raw.csv"
    params:
        config["metad_color"]        
    output:
        "results/cov_norm_good.csv",
        "results/cov_global_good.svg",
        "results/cov_median_good.svg",
        "results/cov_mean_good.svg",
        "results/cov_norm_raw.csv",
        "results/cov_global_raw.svg",
        "results/cov_median_raw.svg",
        "results/cov_mean_raw.svg",
    conda:
        "envs/r.yaml"
    log:
        "logs/coverage/stats_plot.log"    
    script:
        "scripts/cov_stats_all.R"

rule mapq_distribution:
    input:
        "analysis/{sample}/mapq.csv",
        "files/chromosome_names.csv"
    output:
        "analysis/{sample}/mapq_distribution.svg"
    conda:
        "envs/r.yaml"
    log:
        "logs/mapq-dist/{sample}.log"
    script:
        "scripts/mapq-distribution.R"

rule cov_distribution:
    input:
        "analysis/{sample}/cov.csv",
        "files/chromosome_names.csv"
    output:
        "analysis/{sample}/cov_distribution.svg"
    conda:
        "envs/r.yaml"
    log:
        "logs/cov-dist/{sample}.log"
    script:
        "scripts/coverage-distribution.R"

# rule plot_bamstats:
#     input:
#         "analysis/{sample}/snps.bam.stats"
#     output:
#         directory("analysis/{sample}/bamstats")
#     conda:
#         "envs/plot-bamstats.yaml"
#     log:
#         "logs/plot-bamstats/{sample}.log"
#     shell:
#         "plot-bamstats -p {output}/ {input} &> {log}"

rule mapped_plot:
    input:
        "results/mapping_stats.txt",
        config["sample_file"]
    output:
        "results/mapped_reads.svg"
    conda:
        "envs/r.yaml"
    log:
        "logs/stats/mapped.log"
    script:
        "scripts/mapped_reads.R"

rule mapq_plot:
    input:
        "analysis/{sample}/mapq_window.bed",
        "files/chromosome_names.csv",
        "files/loci_to_plot.tsv"
    output:
        "analysis/{sample}/mapq.svg"
    conda:
        "envs/r.yaml"
    log:
        "logs/mapq_plot/{sample}.log"
    script:
        "scripts/mapq.R"