
# =================================================================================================
#   Per sample | Run Mosdepth to get depth per window of all (raw) reads
# =================================================================================================

rule mosdepth:
    input:
        bam = SAMPLES_DIR / "snippy" / "{sample}" / "snps.bam",
        bai = SAMPLES_DIR / "snippy" / "{sample}" / "snps.bam.bai"
    output:
        bed = INT_SAMPLES_DIR / "mosdepth" / "{sample}" / "coverage.regions.bed.gz"
    params:
        window = config["depth_quality"]["mosdepth"]["window"],
        extra = config["depth_quality"]["mosdepth"]["extra"],
        outdir = INT_SAMPLES_DIR / "mosdepth"
    conda: 
        "../envs/depth.yaml"
    threads:
      config["depth_quality"]["mosdepth"]["threads"]    
    log:
        "logs/samples/mosdepth/mosdepth_{sample}.log"
    shell:
        "mosdepth -n --by {params.window} -t {threads} {params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage {input.bam} "
        "&> {log}"

# =================================================================================================
#   Per sample | Get mean MAPQ by window
# =================================================================================================

rule mapq:
    input:
       SAMPLES_DIR / "snippy" / "{sample}" / "snps.bam",
       rules.mosdepth.output.bed
    output:
        bed = INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq.bed",
        window_bed = INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq_window.bed"
    conda:
        "../envs/samtools.yaml"
    log:
       "logs/samples/depth_quality/mapq_{sample}.log"
    script:
        "../scripts/pileup_mapq.sh"
# =================================================================================================
#   Per sample | Intersect MAPQ with depth by window and GFF
# =================================================================================================
rule mapq_depth:
    input:
        mapqbed = rules.mapq.output.window_bed,
        depthbed = rules.mosdepth.output.bed,
        gff = rules.move_liftoff.output,
        global_mode = SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv"
    output:
        depthmapq = SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq_depth_window.bed",
        tsv = SAMPLES_DIR / "depth_quality" / "{sample}" / "feature_mapq_depth.tsv"
    conda:
        "../envs/samtools.yaml"
    log: 
        "logs/samples/depth_quality/mapq_depth_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapq_depth.xsh "
        "-mi {input.mapqbed} "
        "-di {input.depthbed} "
        "-gi {input.gff} "
        "-gmi {input.global_mode} "
        "-sp {wildcards.sample} "
        "-dmo {output.depthmapq} "
        "-o {output.tsv} &> {log}"
