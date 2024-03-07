# Get the coverage of all the mapped reads per window along all chromosomes
rule mosdepth:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai
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
        "{params.outdir}/{wildcards.sample}/coverage {input.bam} "
        "&> {log}"

# Get the coverage of the good quality (above a MAPQ value) mapped reads per window along all chromosomes
rule mosdepth_good:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai
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
        "logs/mosdepth/good_{sample}.log"
    shell:
        "mosdepth -n --by {params.window} --mapq {params.min_mapq} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage_good {input.bam} "
        "&> {log}"

# Get bam files with only good alignments
rule bam_good:
    input:
        bam = rules.snippy.output.bam
    output:
        bam_good = OUTDIR / "samtools" / "{sample}" / "snps_good.bam",
        bai_good = OUTDIR / "samtools" / "{sample}" / "snps_good.bam.bai"
    conda:
        "../envs/samtools.yaml"
    params:
        min_mapq = config["coverage_quality"]["mosdepth"]["min_mapq"]   
    log:
        "logs/stats/bam_good_{sample}.log"
    shell:
        "samtools view -q {params.min_mapq} -b {input} > {output.bam_good} 2> {log} && "
        "samtools index {output.bam_good} -o {output.bai_good} 2>> {log} "

# Get distribution of MAPQ and Coverage values in all alignments and only good alignments
rule samtools_stats:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai,
        bam_good = rules.bam_good.output.bam_good,
        bai_good = rules.bam_good.output.bai_good,
        ref = rules.snippy.output.ref
    output:
        mapq = OUTDIR / "samtools" / "{sample}" / "distrib_mapq.csv",
        cov = OUTDIR / "samtools" / "{sample}" / "distrib_cov.csv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/stats/samtools_stats_{sample}.log"
    shell:
        "xonsh workflow/scripts/samtools-stats.xsh {wildcards.sample} {input.bam} {input.bam_good} {input.ref} {output.mapq} {output.cov} &> {log}"

# Run samtools stats on BAM with all alignments
rule bamstats:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai
    output:
        stats = OUTDIR / "samtools" / "{sample}" / "snps.bam.stats",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/stats/bamstats_{sample}.log"
    shell:
        "samtools stats {input.bam} 1> {output.stats} 2> {log}"
# Get stats on number of mapped reads
rule mapped_edit:
    input:
        stats = rules.bamstats.output.stats 
    output: 
        mapstats = OUTDIR / "samtools" / "{sample}" / "mapping_stats.txt"
    log:
        "logs/stats/mapped_edit_{sample}.log"
    shell:
        "grep reads {input.stats} | cut -d'#' -f1 | cut -f 2- | grep . > {output.mapstats} 2> {log} "
        " && "
        'sed -i "s/$/:\\{wildcards.sample}/" {output.mapstats} 2>> {log}'

# Join the mapping stats of all samples
rule mapped_cat:
    input:
        expand(rules.mapped_edit.output.mapstats, sample=SAMPLES)   
    output: 
        stats = DATASET_OUTDIR / "files" / "mapping_stats.txt"
    log:
        "logs/stats/mapped_cat.log"
    shell:
       'cat {input} > {output}'  

# Get the MAPQ per position and per window
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

# Add the MAPQ and Coverage to the gff file
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
        "logs/mapq/mapqcov2gff_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapqcov2gff.xsh {input.mapqbed} {input.covbed} {input.gff} {output.covmapq} {output.newgff} &> {log}"

# Run RepeatModeler for each reference genome
rule repeat_modeler:
    input:
        rules.links.output
    output:
        known = REFDIR / "{lineage}" / "repeats" / "{lineage}_known.fa",
        unknown = REFDIR / "{lineage}" / "repeats" / "{lineage}_unknown.fa"
    params:
        repdir = "repeats"
    threads:
        config["coverage_quality"]["repeats"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/repeats/repeatmodeler_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-modeler.sh {threads} {input} {params.repdir} &> {log}"

# Run RepeatMasker for each reference genome. Obtain a BED file with the location of the reapeat sequences
rule repeats:
    input:
        database = config["coverage_quality"]["repeats"]["repeats_database"],
        fasta = rules.links.output,
        known = rules.repeat_modeler.output.known,
        unknown = rules.repeat_modeler.output.unknown
    output:
        REFDIR / "{lineage}" / "repeats" / "{lineage}_repeats.bed"
    threads:
        config["coverage_quality"]["repeats"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/repeats/repeatmasker_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-masker.sh {threads} {input.database} {input.fasta} {input.known} {input.unknown} {output} &> {log}"

# Get coverage stats for each window and each chromosome
def raw_coverage_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "coverage": OUTDIR / "mosdepth" / s["sample"] / "coverage.regions.bed.gz" ,
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
    }
rule raw_coverage:
    input:
        unpack(raw_coverage_input)
    output:
        chromosome = OUTDIR / "mosdepth" / "{sample}" / "raw_chromosome_coverage.tsv",
        regions = OUTDIR / "mosdepth" / "{sample}" / "raw_regions_coverage.tsv",
        structure = OUTDIR / "mosdepth" / "{sample}" / "raw_structural_variants.tsv"
    params:
        region = config["coverage_quality"]["mosdepth"]["window"],
        smooth = config["coverage_quality"]["ploidy"]["smoothing_size"],
        repeats_threshold = config["coverage_quality"]["repeats"]["repeats_fraction"],
        change_threshold = config["coverage_quality"]["ploidy"]["change"],
        repeat_category_threshold = config["coverage_quality"]["repeats"]["category_fraction"]
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/coverage/{sample}.log"
    shell:
        "xonsh workflow/scripts/coverage_analysis.xsh "
        "-b {input.coverage} "
        "-rp {input.repeats} "
        "-ch {output.chromosome} "
        "-rg {output.regions} "
        "-sv {output.structure} "
        "-sn {wildcards.sample} "
        "-rs {params.region} "
        "-ss {params.smooth} "
        "-rt {params.repeats_threshold} "
        "-ct {params.change_threshold} "
        "-rct {params.repeat_category_threshold} "
        "&> {log}"

# Get coverage stats for each window and each chromosome for the good quality mappings
def good_coverage_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "coverage": OUTDIR / "mosdepth" / s["sample"] / "coverage_good.regions.bed.gz" ,
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
    }
rule good_coverage:
    input:
        unpack(good_coverage_input)
    output:
        chromosome = OUTDIR / "mosdepth" / "{sample}" / "good_chromosome_coverage.tsv",
        regions = OUTDIR / "mosdepth" / "{sample}" / "good_regions_coverage.tsv",
        structure = OUTDIR / "mosdepth" / "{sample}" / "good_structural_variants.tsv"
    params:
        region = config["coverage_quality"]["mosdepth"]["window"],
        smooth = config["coverage_quality"]["ploidy"]["smoothing_size"],
        repeats_threshold = config["coverage_quality"]["repeats"]["repeats_fraction"],
        change_threshold = config["coverage_quality"]["ploidy"]["change"],
        repeat_category_threshold = config["coverage_quality"]["repeats"]["category_fraction"]
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/coverage/{sample}.log"
    shell:
        "xonsh workflow/scripts/coverage_analysis.xsh "
        "-b {input.coverage} "
        "-rp {input.repeats} "
        "-ch {output.chromosome} "
        "-rg {output.regions} "
        "-sv {output.structure} "
        "-sn {wildcards.sample} "
        "-rs {params.region} "
        "-ss {params.smooth} "
        "-rt {params.repeats_threshold} "
        "-ct {params.change_threshold} "
        "-rct {params.repeat_category_threshold} "
        "&> {log}"
        
# Get the coverage stats of all samples
rule dataset_coverage:
    input:
        g = expand(rules.good_coverage.output.chromosome,sample=SAMPLES),
        r = expand(rules.raw_coverage.output.chromosome,sample=SAMPLES),
        sv = expand(rules.good_coverage.output.structure,sample=SAMPLES)
    output:
        allg = DATASET_OUTDIR / "files" / "coverage_good.tsv",
        allr = DATASET_OUTDIR / "files" / "coverage_raw.tsv",
        allsv = DATASET_OUTDIR / "files" / "structural_variants.tsv"
    log:
        "logs/coverage/dataset_coverage.log"
    script:
        "../scripts/dataset_coverage.sh"