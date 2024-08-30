
# =================================================================================================
#   Per sample | Run Mosdepth to get depth per window of all (raw) reads
# =================================================================================================

rule mosdepth:
    input:
        bam = OUTDIR / "snippy" / "{sample}" / "snps.bam",
        bai = OUTDIR / "snippy" / "{sample}" / "snps.bam.bai"
    output:
        bed = OUTDIR / "mosdepth" / "{sample}" / "coverage.regions.bed.gz"
    params:
        window = config["depth_quality"]["mosdepth"]["window"],
        extra = config["depth_quality"]["mosdepth"]["extra"],
        outdir = OUTDIR / "mosdepth"
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
       OUTDIR / "snippy" / "{sample}" / "snps.bam",
       rules.mosdepth.output.bed
    output:
        bed = OUTDIR / "depth_quality" / "{sample}" / "mapq.bed",
        window_bed = OUTDIR / "depth_quality" / "{sample}" / "mapq_window.bed"
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
        gff = rules.liftoff.output.polished,
        global_mode = OUTDIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv"
    output:
        depthmapq = OUTDIR / "depth_quality" / "{sample}" / "mapq_depth_window.bed",
        tsv = OUTDIR / "depth_quality" / "{sample}" / "feature_mapq_depth.tsv"
    conda:
        "../envs/samtools.yaml"
    log: 
        "logs/samples/depth_quality/mapq_depth_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapq_depth.xsh -m {input.mapqbed} -c {input.depthbed} -g {input.gff} -cm {output.depthmapq} -gm {input.global_mode} -s {wildcards.sample} -o {output.tsv} &> {log}"
