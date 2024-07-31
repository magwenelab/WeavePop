# =================================================================================================
#   Per sample | Intersect MAPQ with depth by region and GFF
# =================================================================================================
rule mapq_depth:
    input:
        mapqbed = rules.mapq.output.region_bed,
        depthbed = rules.mosdepth.output.bed,
        gff = rules.liftoff.output.polished,
        mode = rules.depth_distribution.output.global_mode
    output:
        depthmapq = OUTDIR / "depth_quality" / "{sample}" / "mapq_depth_region.bed",
        tsv = OUTDIR / "depth_quality" / "{sample}" / "feature_mapq_depth.tsv"
    conda:
        "../envs/samtools.yaml"
    log: 
        "logs/samples/depth_quality/mapq_depth_{sample}.log"
    shell:
        "xonsh workflow/scripts/mapq_depth.xsh -m {input.mapqbed} -c {input.depthbed} -g {input.gff} -cm {output.depthmapq} -gm {input.mode} -s {wildcards.sample} -o {output.tsv} &> {log}"
