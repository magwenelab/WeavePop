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

# Get 
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

rule mapped_cat:
    input:
        expand(rules.mapped_edit.output.mapstats, sample=SAMPLES)   
    output: 
        stats = DATASET_OUTDIR / "mapping_stats.txt"
    log:
        "logs/stats/mapped_cat.log"
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
        "logs/mapq/mapqcov2gff_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapqcov2gff.xsh {input.mapqbed} {input.covbed} {input.gff} {output.covmapq} {output.newgff} &> {log}"

rule coverage:
    input:
        rules.mosdepth.output.bed,
        rules.mosdepth_good.output.bed,
        CHROM_NAMES,
    output:
        good = OUTDIR / "mosdepth" / "{sample}" / "good_stats_regions.tsv",
        raw = OUTDIR / "mosdepth" / "{sample}" / "raw_stats_regions.tsv"
    conda:
        "../envs/r.yaml"
    log:
        "logs/coverage/{sample}.log"
    script:
        "../scripts/coverage.R"

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
        config["coverage_quality"]["ploidy"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/ploidy/repeatmodeler_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-modeler.sh {threads} {params.lin} {params.refdir}/{params.lin}/{params.repdir} &> {log}"


rule repeats:
    input:
        REFDIR / "{lineage}" / "repeats" / "{lineage}_known.fa",
        REFDIR / "{lineage}" / "repeats" / "{lineage}_unknown.fa"
    output:
        REFDIR / "{lineage}" / "repeats" / "05_full_out" / "{lineage}.full_mask.bed"
    params:
        refdir = REFDIR ,
        lin = "{lineage}",
        repdir = "repeats",
        database = config["coverage_quality"]["ploidy"]["repeats_database"]
    threads:
        config["coverage_quality"]["ploidy"]["repeats_threads"]
    conda:
        "../envs/repeatmasker.yaml"
    log:
        "logs/ploidy/repeatmasker_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-masker.sh {threads} {params.lin} {params.refdir}/{params.lin}/{params.repdir} {params.database} &> {log}"

rule smoothing:
    input:
        rules.coverage.output.good
    output:
        OUTDIR / "mosdepth" / "{sample}" / "smooth_good_stats_regions.tsv"
    params:
        size = config["coverage_quality"]["ploidy"]["smoothing_size"]
    log:
        "logs/ploidy/smoothing_{sample}.log"
    script:
        "../scripts/median_filtering.py"

rule intersect:
    input:
        unpack(intersect_input)
    output:
        OUTDIR / "mosdepth" / "{sample}" / "smooth_repeats_good_stats_regions.tsv"
    conda:
        "../envs/repeatmasker.yaml"
    params:
        dir = OUTDIR / "mosdepth",
        sample = "{sample}",
        file = "intersect.bed"
    log: 
        "logs/ploidy/intersect_{sample}.log"
    shell:
        """
        tail -n +2 {input.sampletsv} | bedtools intersect -wa -c -a stdin -b {input.maskbed} > {params.dir}/{params.sample}/{params.file} 2> {log}
        head -n1 {input.sampletsv} | paste - <(echo "Nb_Repeats") | cat - {params.dir}/{params.sample}/{params.file} > {output} 2> {log}
        """

rule ploidy_table:
    input:
        rules.intersect.output
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