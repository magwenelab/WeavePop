configfile: "config.yaml"

import pandas as pd

samplefile=(pd.read_csv(config["sample_file"], sep=","))
samples=list(set(samplefile["sample"]))
LINS=list(set(samplefile["group"]))
REFDIR = str(config["reference_directory"])


rule all:
    input:
        expand("analysis/{sample}/coverage.regions.bed.gz",sample=samples),
        expand("analysis/{sample}/coverage_good.regions.bed.gz",sample=samples),
        expand("analysis/{sample}/mapq.csv",sample=samples),
        expand("analysis/{sample}/cov.csv",sample=samples),
        expand("analysis/{sample}/snps.bam.stats",sample=samples),
        "results/mapping_stats.txt",
        expand("analysis/{sample}/mapq_window.bed",sample=samples),
        expand("analysis/{sample}/mapq.bed",sample=samples),
        expand("analysis/{sample}/mapq_window.bed",sample=samples),
        expand("analysis/{sample}/mapq_cov_window.bed",sample=samples),
        expand("analysis/{sample}/annotation.gff",sample=samples),
        expand(REFDIR + "{lineage}.gff.tsv", lineage= LINS),
        config["locitsv"]

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

rule samtools_stats:
    input:
        bam = "analysis/{sample}/snps.bam",
        ref = "analysis/{sample}/ref.fa"
    output:
        mapq = "analysis/{sample}/mapq.csv",
        cov = "analysis/{sample}/cov.csv"
    conda: 
        "envs/samtools.yaml"
    log:
        "logs/stats/{sample}.log"
    shell:
        "xonsh scripts/samtools-stats.xsh {wildcards.sample} {input.bam} {input.ref} {output.mapq} {output.cov} &> {log}"

rule bamstats:
    input:
        "analysis/{sample}/snps.bam"
    output:
        "analysis/{sample}/snps.bam.stats"
    conda:
        "envs/samtools.yaml"
    log:
        "logs/bamstats/{sample}.log"
    shell:
        "samtools stats {input} 1> {output} 2> {log}"

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

rule mapq:
    input:
       "analysis/{sample}/snps.bam",
       "analysis/{sample}/coverage.regions.bed.gz"
    output:
        "analysis/{sample}/mapq.bed",
        "analysis/{sample}/mapq_window.bed" 
    conda:
        "envs/samtools.yaml"
    log:
       "logs/mapq/{sample}.log"
    script:
        "scripts/pileup_mapq.sh"

rule mapqcov2gff:
    input:
        mapqbed = "analysis/{sample}/mapq_window.bed",
        covbed = "analysis/{sample}/coverage.regions.bed.gz",
        gff = "analysis/{sample}/lifted.gff_polished"
    output:
        covmapq = "analysis/{sample}/mapq_cov_window.bed",
        newgff = "analysis/{sample}/annotation.gff"
    conda:
        "envs/samtools.yaml"
    log: 
        "logs/gff/{sample}.log"
    shell:
        "xonsh scripts/mapqcov2gff.xsh {input.mapqbed} {input.covbed} {input.gff} {output.covmapq} {output.newgff} &> {log}"

rule gff2tsv:
    input:
        REFDIR + "{lineage}.gff"
    output:
        REFDIR + "{lineage}.gff.tsv"
    conda:
        "envs/agat.yaml"
    log:
        "logs/references/{lineage}_gff2tsv.log"
    shell:
        "agat_convert_sp_gff2tsv.pl -gff {input} -o {output} "
        "&> {log} && "
        "rm {wildcards.lineage}.agat.log || true"

rule loci:
    input:
        expand(REFDIR + "{lineage}.gff.tsv", lineage=LINS)
    output:
        config["locitsv"]
    params:
        loci=config["loci"]
    log: 
        "logs/references/loci.log"
    run:
        if config["loci"] == "":
            shell("touch {output}")
        else:
            shell("xonsh scripts/loci.xsh {params.loci} -o {output} {input} &> {log}")