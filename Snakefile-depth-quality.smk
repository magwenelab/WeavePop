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
        expand("analysis/{sample}/bamstats", sample=samples),
        "results/mapped_reads.svg",
        expand("analysis/{sample}/mapq_window.bed", sample=samples),
        expand("analysis/{sample}/mapq.svg", sample=samples),
        expand("analysis/{sample}/annotation.gff", sample=samples),

rule mosdepth:
    input:
        "analysis/{sample}/snps.bam"
    output:
        "analysis/{sample}/coverage.regions.bed.gz"
    params:
        window = config["mosdepth_window"]
    conda: 
        "envs/depth.yaml"
    threads:
       config["threads_mosdepth"]     
    log:
        "logs/mosdepth/{sample}.log"
    shell:
        "mosdepth -n --by {params.window} -t {threads} "
        "analysis/{wildcards.sample}/coverage {input} "
        "&> {log}"

rule mosdepth_good:
    input:
        "analysis/{sample}/snps.bam"
    output:
        "analysis/{sample}/coverage_good.regions.bed.gz"
    params:
        window = config["mosdepth_window"],
        min_mapq = config["mosdepth_min_mapq"]
    conda: 
        "envs/depth.yaml"
    threads:
       config["threads_mosdepth"]     
    log:
        "logs/mosdepth_good/{sample}.log"
    shell:
        "mosdepth -n --by {params.window} --mapq {params.min_mapq} -t {threads} "
        "analysis/{wildcards.sample}/coverage_good {input} "
        "&> {log}"

rule coverage_plot:
    input:
        "analysis/{sample}/coverage.regions.bed.gz",
        "analysis/{sample}/coverage_good.regions.bed.gz",
        "files/chromosome_names.csv",
        config["locitsv"]
    output:
        "analysis/{sample}/coverage.svg",
        "analysis/{sample}/coverage_stats.svg",
        "analysis/{sample}/coverage_stats_good.csv",
        "analysis/{sample}/coverage_stats_raw.csv",
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
    log:
        "logs/coverage/stats_plot.log"    
    script:
        "scripts/cov_stats_all.R"

rule samtools_stats:
    input:
        bam = "analysis/{sample}/snps.bam",
        ref = "analysis/{sample}/ref.fa"
    output:
        mapq = "analysis/{sample}/mapq.csv",
        cov = "analysis/{sample}/cov.csv"
    log:
        "logs/stats/{sample}.log"
    shell:
        "xonsh scripts/samtools-stats.xsh {wildcards.sample} {input.bam} {input.ref} {output.mapq} {output.cov} &> {log}"

rule mapq_distribution:
    input:
        "analysis/{sample}/mapq.csv",
        "files/chromosome_names.csv"
    output:
        "analysis/{sample}/mapq_distribution.svg"
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
    log:
        "logs/cov-dist/{sample}.log"
    script:
        "scripts/coverage-distribution.R"

rule bamstats:
    input:
        "analysis/{sample}/snps.bam"
    output:
        "analysis/{sample}/snps.bam.stats"
    conda:
        "envs/depth.yaml"
    log:
        "logs/bamstats/{sample}.log"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"

rule plot_bamstats:
    input:
        "analysis/{sample}/snps.bam.stats"
    output:
        directory("analysis/{sample}/bamstats")
    conda:
        "envs/depth.yaml"
    log:
        "logs/plot-bamstats/{sample}.log"
    shell:
        "plot-bamstats -p {output}/ {input} &> {log}"

rule mapped_edit:
    input:
        "analysis/{sample}/snps.bam.stats" 
    output: 
        temp("analysis/{sample}/mapping_stats.txt")
    shell:
        "grep reads {input} | cut -d'#' -f1 | cut -f 2- | grep . > {output} "
        " && "
        'sed -i "s/$/:\\{wildcards.sample}/" {output}'

rule mapped_cat:
    input:
        expand("analysis/{sample}/mapping_stats.txt", sample=samples)   
    output: 
        "results/mapping_stats.txt"
    shell:
       'cat {input} > {output}'  

rule mapped_plot:
    input:
        "results/mapping_stats.txt",
        config["sample_file"]
    output:
        "results/mapped_reads.svg"
    log:
        "logs/stats/mapped.log"
    script:
        "scripts/mapped_reads.R"

rule mapq:
   input:
       "analysis/{sample}/snps.bam",
       "analysis/{sample}/coverage.regions.bed.gz"
   output:
        "analysis/{sample}/mapq.bed",
        "analysis/{sample}/mapq_window.bed" 
   log:
       "logs/mapq/{sample}.log"
   script:
        "scripts/pileup_mapq.sh"

rule mapq_plot:
    input:
        "analysis/{sample}/mapq_window.bed",
        "files/chromosome_names.csv",
        config["locitsv"]
    output:
        "analysis/{sample}/mapq.svg"
    log:
        "logs/mapq_plot/{sample}.log"
    script:
        "scripts/mapq.R"

rule mapqcov2gff:
    input:
        mapqbed = "analysis/{sample}/mapq_window.bed",
        covbed = "analysis/{sample}/coverage.regions.bed.gz",
        gff = "analysis/{sample}/lifted.gff_polished"
    output:
        covmapq = "analysis/{sample}/mapq_cov_window.bed",
        newgff = "analysis/{sample}/annotation.gff"
    log: 
        "logs/gff/{sample}.log"
    shell:
        "xonsh scripts/mapqcov2gff.xsh {input.mapqbed} {input.covbed} {input.gff} {output.covmapq} {output.newgff} &> {log}"