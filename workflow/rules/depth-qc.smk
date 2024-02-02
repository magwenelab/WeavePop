# Get the coverage of all the mapped reads per window along all chromosomes
rule mosdepth:
    input:
        bam = rules.snippy.output.bam
    output:
        bed = OUTDIR / "mosdepth" / "{sample}" / "coverage.regions.bed.gz"
    params:
        window = config["coverage_quality"]["mosdepth"]["window"],
        extra = config["coverage_quality"]["mosdepth"]["extra"],
        outdir = OUTDIR / "mosdepth"
    conda: 
        "../envs/depth.yaml"
    threads:
      config["coverage_quality"]["mosdepth"]["threads"]    
    log:
        "logs/mosdepth/{sample}.log"
    shell:
        "mosdepth -n --by {params.window} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage {input} "
        "&> {log}"

# Get the coverage of the good quality (above a MAPQ value) mapped reads per window along all chromosomes
rule mosdepth_good:
    input:
        bam = rules.snippy.output.bam
    output:
        bed = OUTDIR / "mosdepth" / "{sample}" / "coverage_good.regions.bed.gz"
    params:
        window = config["coverage_quality"]["mosdepth"]["window"],
        extra = config["coverage_quality"]["mosdepth"]["extra"],
        min_mapq = config["coverage_quality"]["mosdepth"]["min_mapq"],
        outdir = OUTDIR / "mosdepth"
    conda: 
        "../envs/depth.yaml"
    threads:
       config["coverage_quality"]["mosdepth"]["threads"]   
    log:
        "logs/mosdepth_good/{sample}.log"
    shell:
        "mosdepth -n --by {params.window} --mapq {params.min_mapq} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage_good {input} "
        "&> {log}"

# Get 
rule samtools_stats:
    input:
        bam = rules.snippy.output.bam,
        ref = rules.snippy.output.ref
    output:
        mapq = OUTDIR / "samtools" / "{sample}" / "distrib_mapq.csv",
        cov = OUTDIR / "samtools" / "{sample}" / "distrib_cov.csv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/stats/{sample}.log"
    shell:
        "xonsh workflow/scripts/samtools-stats.xsh {wildcards.sample} {input.bam} {input.ref} {output.mapq} {output.cov} &> {log}"

rule bamstats:
    input:
        bam = rules.snippy.output.bam
    output:
        stats = OUTDIR / "samtools" / "{sample}" / "snps.bam.stats",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bamstats/{sample}.log"
    shell:
        "samtools stats {input.bam} 1> {output.stats} 2> {log}"

rule mapped_edit:
    input:
        stats = rules.bamstats.output.stats 
    output: 
        mapstats = OUTDIR / "samtools" / "{sample}" / "mapping_stats.txt"
    shell:
        "grep reads {input.stats} | cut -d'#' -f1 | cut -f 2- | grep . > {output.mapstats} "
        " && "
        'sed -i "s/$/:\\{wildcards.sample}/" {output.mapstats}'

rule mapped_cat:
    input:
        expand(rules.mapped_edit.output.mapstats, sample=SAMPLES)   
    output: 
        stats = DATASET_OUTDIR / "mapping_stats.txt"
    shell:
       'cat {input} > {output}'  

rule mapq:
    input:
       rules.snippy.output.bam,
       rules.mosdepth.output.bed
    output:
        bed = OUTDIR / "samtools" / "{sample}" / "mapq.bed",
        winbed = OUTDIR / "samtools" / "{sample}" / "mapq_window.bed" 
    conda:
        "../envs/samtools.yaml"
    log:
       "logs/mapq/{sample}.log"
    script:
        "../scripts/pileup_mapq.sh"

rule mapqcov2gff:
    input:
        mapqbed = rules.mapq.output.winbed,
        covbed = rules.mosdepth.output.bed,
        gff = rules.liftoff.output.polished
    output:
        covmapq = OUTDIR / "samtools" / "{sample}" / "mapq_cov_window.bed",
        newgff = OUTDIR / "samtools" / "{sample}" / "annotation.gff"
    conda:
        "../envs/samtools.yaml"
    log: 
        "logs/gff/{sample}.log"
    shell:
        "xonsh workflow/scripts/mapqcov2gff.xsh {input.mapqbed} {input.covbed} {input.gff} {output.covmapq} {output.newgff} &> {log}"

rule coverage:
    input:
        rules.mosdepth.output.bed,
        rules.mosdepth_good.output.bed,
        CHROM_NAMES,
    output:
        good = OUTDIR / "mosdepth" / "{sample}" / "good_stats_regions.csv",
        raw = OUTDIR / "mosdepth" / "{sample}" / "raw_stats_regions.csv"
    conda:
        "../envs/r.yaml"
    log:
        "logs/coverage/{sample}.log"
    script:
        "../scripts/coverage.R"

rule ploidy_table:
    input:
        rules.coverage.output.good
    output:
        OUTDIR / "mosdepth" / "{sample}" / "ploidy_table.csv"
    conda:
        "../envs/r.yaml"
    params:
        diff_threshold = config["coverage_quality"]["ploidy"]["diff"],  
        size_threshold = config["coverage_quality"]["ploidy"]["size"], 
        change_threshold = config["coverage_quality"]["ploidy"]["change"]  
    log:
        "logs/ploidy/{sample}.log"
    script:
        "../scripts/ploidy_table.R"