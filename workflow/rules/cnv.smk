# =================================================================================================
#   Per sample | Run Mosdepth to get depth per window of good quality reads
# =================================================================================================


rule mosdepth_good:
    input:
        bam=SAMPLES_DIR / "snippy" / "{sample}" / "snps.bam",
        bai=SAMPLES_DIR / "snippy" / "{sample}" / "snps.bam.bai",
    output:
        bed=INT_SAMPLES_DIR / "mosdepth" / "{sample}" / "coverage_good.regions.bed.gz",
    params:
        window=config["depth_quality"]["mosdepth"]["window"],
        extra=config["depth_quality"]["mosdepth"]["extra"],
        min_mapq=config["depth_quality"]["flag_quality"]["min_mapq"],
        outdir=INT_SAMPLES_DIR / "mosdepth",
    log:
        LOGS / "samples" / "depth_quality" / "mosdepth_good_{sample}.log",
    threads: config["depth_quality"]["mosdepth"]["threads"]
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/depth.yaml"
    shell:
        "mosdepth "
        "-n "
        "--by {params.window} "
        "--mapq {params.min_mapq} "
        "-t {threads} "
        "{params.extra} "
        "{params.outdir}/{wildcards.sample}/coverage_good "
        "{input.bam} "
        "&> {log}"


# =================================================================================================
#   Per sample | Normalize and smooth depth by windows
# =================================================================================================


rule depth_by_windows:
    input:
        depth=rules.mosdepth_good.output.bed,
        genome_wide_depth=SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_chrom_good.tsv",
    output:
        INT_SAMPLES_DIR / "depth_quality" / "{sample}" / "depth_by_windows.tsv",
    params:
        smoothing_size=config["cnv"]["smoothing_size"],
    log:
        LOGS / "samples" / "depth_quality" / "depth_by_windows_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/depth_by_windows.xsh "
        "-di {input.depth} "
        "-gi {input.genome_wide_depth} "
        "-do {output} "
        "-s {params.smoothing_size} "
        "&> {log}"


# =================================================================================================
#   Per sample | Call CNVs
# =================================================================================================


rule cnv_calling:
    input:
        unpack(cnv_calling_input),
    output:
        SAMPLES_DIR / "cnv" / "{sample}" / "cnv_calls.tsv",
    params:
        window_size=config["depth_quality"]["mosdepth"]["window"],
        depth_threshold=config["cnv"]["depth_threshold"],
    log:
        LOGS / "samples" / "cnv" / "cnv_calling_{sample}.log",
    resources:
        tmpdir=TEMPDIR,
    conda:
        "../envs/samtools.yaml"
    shell:
        "xonsh workflow/scripts/cnv_calling.xsh "
        "-di {input.depth} "
        "-ri {input.repeats} "
        "-co {output} "
        "-sp {wildcards.sample} "
        "-wp {params.window_size} "
        "-dp {params.depth_threshold} &> {log}"
