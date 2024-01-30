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
    params:
        script = workflow.source_path("../scripts/loci.xsh")
    log: 
        "logs/references/loci.log"
    shell:
        "xonsh {params.script} -g {input.loci} -o {output} {input.refs} &> {log}"

rule coverage_plot:
    input:
        rules.mosdepth.output.bed,
        rules.mosdepth_good.output.bed,
        CHROM_NAMES,
        rules.loci.output.locitable
    output:
        cov = OUTDIR / "plots" / "{sample}" / "coverage.svg",
        stats = OUTDIR / "plots" / "{sample}" / "coverage_stats.svg",
        good = OUTDIR / "files" / "{sample}" / "coverage_stats_good.csv",
        raw = OUTDIR / "files" / "{sample}" / "coverage_stats_raw.csv"
    conda:
        "../envs/r.yaml"
    params:
        script = workflow.source_path("../scripts/coverage.R")
    log:
        "logs/coverage/{sample}.log"
    script:
        "{params.script}"

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
        DATASET_OUTDIR / "plots" / "cov_global_good.svg",
        DATASET_OUTDIR / "plots" / "cov_median_good.svg",
        DATASET_OUTDIR / "plots" / "cov_mean_good.svg",
        DATASET_OUTDIR / "files" / "cov_norm_raw.csv",
        DATASET_OUTDIR / "plots" / "cov_global_raw.svg",
        DATASET_OUTDIR / "plots" / "cov_median_raw.svg",
        DATASET_OUTDIR / "plots" / "cov_mean_raw.svg",
    conda:
        "../envs/r.yaml"
    params:
        config["plotting"]["metadata2color"],
        script = workflow.source_path("../scripts/cov_stats_all.R")
    log:
        "logs/coverage/stats_plot.log"    
    script:
        "{params.script}"

rule mapq_distribution:
    input:
        rules.samtools_stats.output.mapq,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "mapq_distribution.svg"
    conda:
        "../envs/r.yaml"
    params:
        script = workflow.source_path("../scripts/mapq-distribution.R")    
    log:
        "logs/mapq-dist/{sample}.log"
    script:
        "{params.script}"

rule cov_distribution:
    input:
        rules.samtools_stats.output.cov,
        CHROM_NAMES
    output:
        OUTDIR / "plots" / "{sample}" / "cov_distribution.svg"
    conda:
        "../envs/r.yaml"
    params:
        script = workflow.source_path("../scripts/coverage-distribution.R") 
    log:
        "logs/cov-dist/{sample}.log"
    script:
        "{params.script}"

rule mapped_plot:
    input:
        rules.mapped_cat.output.stats,
        SAMPLEFILE
    output:
        DATASET_OUTDIR / "plots" / "mapped_reads.svg"
    conda:
        "../envs/r.yaml"
    log:
        "logs/stats/mapped.log"
    params:
        script = workflow.source_path("../scripts/mapped_reads.R") 
    script:
        "{params.script}"

rule mapq_plot:
    input:
        rules.mapq.output.winbed,
        CHROM_NAMES,
        rules.loci.output.locitable
    output:
        OUTDIR / "plots" / "{sample}" / "mapq.svg"
    conda:
        "../envs/r.yaml"
    params:
        script = workflow.source_path("../scripts/mapq.R")
    log:
        "logs/mapq_plot/{sample}.log"
    script:
        "{params.script}"

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
                
# rule unmapped_count_plot:
#     input:
#         REF_GFF + ".tsv",
#         SAMPLEFILE,
#         expand(REFDIR + "{lineage}_unmapped_features.txt", lineage=LINEAGES)        
#     output:
#         REFDIR + "references_unmapped_count.tsv",
#         REFDIR + "references_unmapped.svg"
#     conda:
#         "../envs/r.yaml"
#     params:
#         script = workflow.source_path("../scripts/count_reference_unmapped.R")
#     log:
#         "logs/references/unmapped_count_plot.log"
#     script:
#         "{params.script}"

# rule unmapped_samples_plot:
#     input:
#         REF_GFF + ".tsv",
#         SAMPLEFILE,
#         expand("analysis/{sample}/unmapped_features.txt", sample=SAMPLES)        
#     output:
#         "results/unmapped_count.tsv",
#         "results/unmapped.svg"
#     conda:
#         "../envs/r.yaml"
#     params:
#         script = workflow.source_path("../scripts/count_sample_unmapped.R")
#     log:
#         "logs/liftoff/unmapped_count_plot.log"
#     script:
#         "{params.script}"