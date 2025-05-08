# =================================================================================================
#   Per sample | Get mean MAPQ by window
# =================================================================================================


rule mapq:
    input:
        INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "snps_good.bam",
        rules.mosdepth_good.output.bed,
    output:
        bed=INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq.bed",
        window_bed=INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq_by_window.bed",
    log:
        LOGS / "samples" / "depth_quality" / "mapq_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/mapq.sh"


# =================================================================================================
#   Per sample | Intersect MAPQ to depth by window and annotated features
# =================================================================================================


rule mapq_depth:
    input:
        mapqbed=rules.mapq.output.window_bed,
        depthbed=rules.mosdepth_good.output.bed,
        gff=rules.reformat_annotation.output.gff,
    output:
        depthmapq=SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq_depth_by_window.bed",
        tsv=SAMPLES_DIR / "depth_quality" / "{sample}" / "mapq_depth_by_feature.tsv",
    log:
        LOGS / "samples" / "depth_quality" / "mapq_depth_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/mapq_depth.xsh "
        "-mi {input.mapqbed} "
        "-di {input.depthbed} "
        "-gi {input.gff} "
        "-sp {wildcards.sample} "
        "-dmo {output.depthmapq} "
        "-o {output.tsv} &> {log}"
