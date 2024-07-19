# =================================================================================================
#   Per sample | Run Mosdepth to get depth per region of all (raw) and good quality reads
# =================================================================================================

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
        "logs/samples/mosdepth/mosdepth_{sample}.log"
    shell:
        "mosdepth -n --by {params.window} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage {input.bam} "
        "&> {log}"

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
        "logs/samples/mosdepth/mosdepth_good_{sample}.log"
    shell:
        "mosdepth -n --by {params.window} --mapq {params.min_mapq} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage_good {input.bam} "
        "&> {log}"

# =================================================================================================
#   Per sample | Get depth by chromosome raw and good
# =================================================================================================

rule depth_by_chrom_raw:
    input:
        rules.mosdepth.output.bed
    output:
        OUTDIR / "mosdepth" / "{sample}" / "depth_by_chrom_raw.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/mosdepth/depth_by_chrom_raw_{sample}.log"
    shell:
        "xonsh workflow/scripts/depth_by_chrom.xsh -b {input} -o {output} -s {wildcards.sample} &> {log}"

rule depth_by_chrom_good:
    input:
        rules.mosdepth_good.output.bed
    output:
        OUTDIR / "mosdepth" / "{sample}" / "depth_by_chrom_good.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/mosdepth/depth_by_chrom_good_{sample}.log"
    shell:
        "xonsh workflow/scripts/depth_by_chrom.xsh -b {input} -o {output} -s {wildcards.sample} &> {log}"

# =================================================================================================
#   Per sample | Get distribution and global mode fo depth (genome-wide depth to normalize)
# =================================================================================================

rule bam_good:
    input:
        bam = rules.snippy.output.bam
    output:
        bam_good = temp(OUTDIR / "samtools" / "{sample}" / "snps_good.bam"),
        bai_good = temp(OUTDIR / "samtools" / "{sample}" / "snps_good.bam.bai")
    conda:
        "../envs/samtools.yaml"
    params:
        min_mapq = config["coverage_quality"]["mosdepth"]["min_mapq"]   
    log:
        "logs/samples/samtools/bam_good_{sample}.log"
    shell:
        "samtools view -q {params.min_mapq} -b {input} > {output.bam_good} 2> {log} && "
        "samtools index {output.bam_good} -o {output.bai_good} 2>> {log} "

def depth_distribution_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "bam": OUTDIR / "snippy" / s["sample"] / "snps.bam" ,
        "bai": OUTDIR / "snippy" / s["sample"] / "snps.bam.bai",
        "bam_good": OUTDIR / "samtools" / s["sample"] / "snps_good.bam",
        "bai_good": OUTDIR / "samtools" / s["sample"] / "snps_good.bam.bai"
        }
rule depth_distribution:
    input:
        unpack(depth_distribution_input)
    output:
        distrib = OUTDIR / "samtools" / "{sample}" / "depth_distribution.tsv",
        global_mode = OUTDIR / "samtools" / "{sample}" / "global_mode.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/samtools/samtools_stats_{sample}.log"
    shell:
        "xonsh workflow/scripts/depth_distribution.xsh -s {wildcards.sample} -b {input.bam} -g {input.bam_good} -do {output.distrib} -go {output.global_mode} &> {log}"

# =================================================================================================
#   Per sample | Normalize depth by chromosome and by region
# =================================================================================================

rule depth_by_chrom_normalized:
    input:
        depth = rules.depth_by_chrom_good.output,
        global_mode = rules.depth_distribution.output.global_mode
    output:
        OUTDIR / "mosdepth" / "{sample}" / "depth_by_chrom_good_normalized.tsv"
    conda: 
        "../envs/samtools.yaml"
    log:
        "logs/samples/mosdepth/depth_by_chrom_good_normalized_{sample}.log"
    shell:
        "xonsh workflow/scripts/depth_by_chrom_normalized.xsh -d {input.depth} -g {input.global_mode} -o {output} &> {log}"

rule depth_by_regions:
    input:
        depth = rules.mosdepth_good.output.bed,
        global_mode = rules.depth_distribution.output.global_mode
    output:
        OUTDIR / "mosdepth" / "{sample}" / "depth_by_regions.tsv"
    conda:
        "../envs/samtools.yaml"
    params:
        smoothing_size = config["coverage_quality"]["cnv"]["smoothing_size"]
    log:
        "logs/samples/mosdepth/depth_by_regions_{sample}.log"
    shell:
        "xonsh workflow/scripts/depth_by_regions.xsh -di {input.depth} -gi {input.global_mode} -do {output} -s {params.smoothing_size} &> {log}"

# =================================================================================================
#   Per lineage | Run RepeatModeler and RepeatMasker
# =================================================================================================

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
        "logs/references/repeats/repeatmodeler_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-modeler.sh {threads} {input} {params.repdir} &> {log}"

rule repeat_masker:
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
        "logs/references/repeats/repeatmasker_{lineage}.log"
    shell:
        "bash workflow/scripts/repeat-masker.sh {threads} {input.database} {input.fasta} {input.known} {input.unknown} {output} &> {log}"

# =================================================================================================
#   Per sample | Intercept depth by regions with repeats and call CNVs
# =================================================================================================

def cnv_calling_input(wildcards):
    s = SAMPLE_REFERENCE.loc[wildcards.sample,]
    return {
        "depth": OUTDIR / "mosdepth" / s["sample"] / "depth_by_regions.tsv",
        "repeats": REFDIR / s["group"]  / "repeats" / (s["group"] + "_repeats.bed")
        }
rule cnv_calling:
    input:
        unpack(cnv_calling_input)
    output:
        OUTDIR / "cnv" / "{sample}" / "cnv_calls.tsv"
    conda:
        "../envs/samtools.yaml"
    params:
        region_size = config["coverage_quality"]["mosdepth"]["window"],
        depth_threshold = config["coverage_quality"]["cnv"]["depth_threshold"]
    log:
        "logs/samples/cnv/cnv_calling_{sample}.log"
    shell:
        "xonsh workflow/scripts/cnv_calling.xsh -di {input.depth} -ri {input.repeats} -co {output} -sp {wildcards.sample} -rp {params.region_size} -dp {params.depth_threshold} &> {log}"

# =================================================================================================
#   Per sample | Get mapping stats
# =================================================================================================

rule mapping_stats:
    input:
        bam = rules.snippy.output.bam,
        bai = rules.snippy.output.bai
    output:
        OUTDIR / "samtools" / "{sample}" / "mapping_stats.tsv"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samples/samtools/mapping_stats_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapping_stats.xsh -b {input.bam} -s {wildcards.sample} -o {output} &> {log}"

# =================================================================================================
#   Per sample | Get mean MAPQ by region
# =================================================================================================

rule mapq:
    input:
       rules.snippy.output.bam,
       rules.mosdepth.output.bed
    output:
        bed = temp(OUTDIR / "samtools" / "{sample}" / "mapq.bed"),
        winbed = OUTDIR / "samtools" / "{sample}" / "mapq_window.bed" 
    conda:
        "../envs/samtools.yaml"
    log:
       "logs/samples/samtools/mapq_{sample}.log"
    script:
        "../scripts/pileup_mapq.sh"

# =================================================================================================
#   Per sample | Intersect MAPQ with depth by region and GFF
# =================================================================================================
rule mapq_depth:
    input:
        mapqbed = rules.mapq.output.winbed,
        depthbed = rules.mosdepth.output.bed,
        gff = rules.liftoff.output.polished,
        mode = rules.depth_distribution.output.global_mode
    output:
        depthmapq = OUTDIR / "samtools" / "{sample}" / "mapq_depth_window.bed",
        tsv = OUTDIR / "samtools" / "{sample}" / "feature_mapq_depth.tsv"
    conda:
        "../envs/samtools.yaml"
    log: 
        "logs/samples/samtools/mapq_depth_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapq_depth.xsh -m {input.mapqbed} -c {input.depthbed} -g {input.gff} -cm {output.depthmapq} -gm {input.mode} -s {wildcards.sample} -o {output.tsv} &> {log}"
