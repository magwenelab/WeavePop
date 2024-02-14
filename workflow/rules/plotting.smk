rule gff2tsv:
    input:
        REFDIR / "{lineage}" / "{lineage}.gff"
    output:
        REFDIR / "{lineage}" / "{lineage}.gff.tsv"
    conda:
        "../envs/agat.yaml"
    log:
        "logs/references/{lineage}_gff2tsv.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} "
        "&> {log} && "
        "rm {wildcards.lineage}.agat.log || true"

rule loci:
    input:
        refs = expand(rules.gff2tsv.output, lineage=LINEAGES),
        loci=LOCI_FILE
    output:
        locitable = DATASET_OUTDIR / "loci_to_plot.tsv"
    log: 
        "logs/references/loci.log"
    shell:
        "xonsh workflow/scripts/loci.xsh -g {input.loci} -o {output} {input.refs} &> {log}"

rule coverage_plot:
    input:
        rules.coverage.output.good,
        rules.coverage.output.raw,
        rules.loci.output.locitable,
        rules.ploidy_table.output
    output:
        cov = OUTDIR / "plots" / "{sample}" / "coverage.png",
        stats = OUTDIR / "plots" / "{sample}" / "coverage_stats.png",
        good = OUTDIR / "mosdepth" / "{sample}" / "good_stats_chroms.csv",
        raw = OUTDIR / "mosdepth" / "{sample}" / "raw_stats_chroms.csv"
    conda:
        "../envs/r.yaml"
    log:
        "logs/coverage_plot/{sample}.log"
    script:
        "../scripts/coverage_plots.R"

rule cat_stats:
    input:
        r = expand(rules.coverage_plot.output.raw,sample=SAMPLES),
        g = expand(rules.coverage_plot.output.good,sample=SAMPLES),
    output:
        allr = DATASET_OUTDIR / "files" / "coverage_raw.csv",
        allg = DATASET_OUTDIR / "files" / "coverage_good.csv",
    log:
        "logs/coverage/cat_stats.log"
    shell:
        "cat {input.r} | grep -v Sample > {output.allr} "
        "&& "
        "cat {input.g} | grep -v Sample > {output.allg} "
        "2> {log}"

rule coverage_stats_plots:
    input:
        SAMPLEFILE,
        rules.cat_stats.output.allg,
        rules.cat_stats.output.allr
    output:
        DATASET_OUTDIR / "files" / "cov_norm_good.csv",
        DATASET_OUTDIR / "plots" / "cov_global_good.png",
        DATASET_OUTDIR / "plots" / "cov_median_good.png",
        DATASET_OUTDIR / "plots" / "cov_mean_good.png",
        DATASET_OUTDIR / "files" / "cov_norm_raw.csv",
        DATASET_OUTDIR / "plots" / "cov_global_raw.png",
        DATASET_OUTDIR / "plots" / "cov_median_raw.png",
        DATASET_OUTDIR / "plots" / "cov_mean_raw.png",
    conda:
        "../envs/r.yaml"
    params:
        config["plotting"]["metadata2color"]
    log:
        "logs/coverage/stats_plot.log"    
    script:
        "../scripts/cov_stats_all.R"

rule mapq_distribution:
    input:
        rules.samtools_stats.output.mapq,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "mapq_distribution.png"
    conda:
        "../envs/r.yaml"   
    log:
        "logs/mapq-dist/{sample}.log"
    script:
        "../scripts/mapq-distribution.R"

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
        "logs/cov-dist/{sample}.log"
    script:
        "../scripts/coverage-distribution.R"

rule mapped_plot:
    input:
        rules.mapped_cat.output.stats,
        SAMPLEFILE
    output:
        DATASET_OUTDIR / "plots" / "mapped_reads.png"
    conda:
        "../envs/r.yaml"
    log:
        "logs/stats/mapped.log"
    script:
        "../scripts/mapped_reads.R"

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
        "logs/mapq_plot/{sample}.log"
    script:
        "../scripts/mapq.R"

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
                