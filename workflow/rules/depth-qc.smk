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
        stats = DATASET_OUTDIR / "mapping_stats.txt"
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

# Get coverage stats for each window and each chromosome
rule coverage:
    input:
        rules.mosdepth.output.bed,
        rules.mosdepth_good.output.bed,
    output:
        good = OUTDIR / "mosdepth" / "{sample}" / "good_stats_regions.tsv",
        raw = OUTDIR / "mosdepth" / "{sample}" / "raw_stats_regions.tsv",
        good_chrom = OUTDIR / "mosdepth" / "{sample}" / "good_stats_chroms.tsv",
        raw_chrom = OUTDIR / "mosdepth" / "{sample}" / "raw_stats_chroms.tsv"
    conda:
        "../envs/r.yaml"
    log:
        "logs/coverage/{sample}.log"
    script:
        "../scripts/coverage.R"

# Get the coverage stats of all samples
rule cat_stats:
    input:
        g = expand(rules.coverage.output.good_chrom,sample=SAMPLES),
        r = expand(rules.coverage.output.raw_chrom,sample=SAMPLES)
    output:
        allg = DATASET_OUTDIR / "files" / "coverage_good.tsv",
        allr = DATASET_OUTDIR / "files" / "coverage_raw.tsv"
    log:
        "logs/coverage/cat_stats.log"
    shell:
        "cat {input.g} | head -n 1 > {output.allg} && "
        "cat {input.r} | head -n 1 > {output.allr} && "
        "tail -q -n +2 {input.g} >> {output.allg} && "
        "tail -q -n +2 {input.r} >> {output.allr}"

# rule repeats_database:
#     output:
#         REFDIR / "repeats_database.fasta"
#     params:
#         dir = REFDIR,
#         database = "https://www.girinst.org/server/RepBase/protected/RepBase29.01.fasta.tar.gz"
#     shell:
#         """
#         wget -O {params.dir}/database {params.database}
#         tar -xvzf {params.dir}/database --one-top-level={params.dir}/database_dir
#         cat {params.dir}/database_dir/*.ref > {params.dir}/database.fasta
#         cat {params.dir}/database.fasta/appendix/*.ref >> {params.dir}/database.fasta
#         rm -rf {params.dir}/database_dir/ {params.dir}/database
#         """

# Run RepeatModeler for each reference genome
rule repeat_modeler:
    input:
        rules.links.output
    output:
        REFDIR / "{lineage}" / "repeats" / "{lineage}_known.fa",
        REFDIR / "{lineage}" / "repeats" / "{lineage}_unknown.fa"
    params:
        refdir = REFDIR ,
        lin = "{lineage}",
        repdir = "repeats"
    threads:
        config["coverage_quality"]["repeats"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/repeats/repeatmodeler_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-modeler.sh {threads} {params.lin} {params.refdir}/{params.lin}/{params.repdir} &> {log}"

# Run RepeatMasker for each reference genome. Obtain a BED file with the location of the reapeat sequences
rule repeats:
    input:
        rules.links.output,
        REFDIR / "{lineage}" / "repeats" / "{lineage}_known.fa",
        REFDIR / "{lineage}" / "repeats" / "{lineage}_unknown.fa"
    output:
        REFDIR / "{lineage}" / "repeats" / "05_full" / "{lineage}.bed"
    params:
        lin = "{lineage}",
        refdir = REFDIR,
        repdir = "repeats",
        database = config["coverage_quality"]["repeats"]["repeats_database"]
    threads:
        config["coverage_quality"]["repeats"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/repeats/repeatmasker_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-masker.sh {threads} {params.lin} {params.refdir}/{params.lin}/{params.repdir} {params.database} &> {log}"

# Get smoothed coverage for each sample
rule smoothing:
    input:
        rules.coverage.output.good
    output:
        OUTDIR / "mosdepth" / "{sample}" / "smooth_coverage_regions.tsv"
    params:
        size = config["coverage_quality"]["ploidy"]["smoothing_size"]
    log:
        "logs/ploidy/smoothing_{sample}.log"
    script:
        "../scripts/median_filtering.py"

# Detect structural variation
rule ploidy_table:
    input:
        rules.smoothing.output
    output:
        OUTDIR / "mosdepth" / "{sample}" / "ploidy_table.tsv"
    conda:
        "../envs/r.yaml"
    params:
        size_threshold = config["coverage_quality"]["ploidy"]["size"], 
        change_threshold = config["coverage_quality"]["ploidy"]["change"]  
    log:
        "logs/ploidy/ploidy_table_{sample}.log"
    script:
        "../scripts/ploidy_table.R"

# Get the fraction of repeated sequences in each window that is a structural variant
rule intersect:
    input:
        unpack(intersect_input)
    output:
        OUTDIR / "mosdepth" / "{sample}" / "structural_variants.tsv"
    conda:
        "../envs/samtools.yaml"
    params:
        threshold = config["coverage_quality"]["repeats"]["repeats_fraction"]
    log: 
        "logs/ploidy/intersect_{sample}.log"
    shell:
        """
        xonsh workflow/scripts/intersect_repeats.xsh -s {input.sampletsv} -r {input.maskbed} -o {output} -t {params.threshold} 2> {log}
        """